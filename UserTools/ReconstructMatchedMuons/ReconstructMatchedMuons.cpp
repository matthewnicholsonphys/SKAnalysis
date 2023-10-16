#include "ReconstructMatchedMuons.h"
#include "fortran_routines.h"
#include "Constants.h"
#include "Algorithms.h"

#include <algorithm>

ReconstructMatchedMuons::ReconstructMatchedMuons():Tool(){}

// This Tool actually does 2 things:
// 1. It writes out relic and muon events that are flagged by the ReconstructMatchedMuons Tool to be written to file
//    (in the process it merges hits from the following AFT if present)
// 2. It performs muon reconstruction for those muons being written out, applying mfmuselect, muboy,
//    and if required, BFF, as well as calculating muon dEdx.
// The reason these two fuctions are merged (ideally muon reconstruction would be an independent Tool)
// is that we don't know until matching is complete whether a muon should be reconstructed and written out.
// This requires re-reading in a past file entry (e.g. we may only know after reading another 100 events
// off disk that the muon has a relic match and should be kept), an operation required both for
// reconstruction and writing the event to file.
// The correct way to do this would be to use a subToolChain that iterates over matched events,
// and has a Tool for reconstruction and a Tool for Writing Out, but since we're lazy we just combine them

bool ReconstructMatchedMuons::Initialise(std::string configfile, DataModel &data){
	
	if(configfile!="")  m_variables.Initialise(configfile);
	//m_variables.Print();
	
	m_data= &data;
	m_log= m_data->Log;
	
	if(!m_variables.Get("verbosity",m_verbose)) m_verbose=1;
	m_variables.Get("noBFF", noBFF);
	
	// get input file reader
	m_variables.Get("rfmReaderName",rfmReaderName);
	rfmReaderLUN = m_data->GetLUN(rfmReaderName);
	if(rfmReaderLUN<0){
		Log(m_unique_name+" Error! Could not find TreeReader '"+rfmReaderName+"' in DataModel"
		    ,v_error,m_verbose);
		return false;
	}
	rfmReader = m_data->Trees.at(rfmReaderName);
	
	// get LUNs for output file writers 
	// (needed to pass common block data from reco algorithms to TTrees etc)
	std::string relicWriterName, muWriterName;
	m_variables.Get("muWriterName", muWriterName);
	m_variables.Get("relicWriterName", relicWriterName);
	muWriterLUN = m_data->GetLUN(muWriterName);
	if(muWriterLUN<0){
		Log(m_unique_name+" Error! Failed to find TreeReader "+muWriterName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	relicWriterLUN = m_data->GetLUN(relicWriterName);
	if(relicWriterLUN<0){
		Log(m_unique_name+" Error! Failed to find "+relicWriterName+" in DataModel!",v_error,m_verbose);
		m_data->vars.Set("StopLoop",1);
		return false;
	}
	
	// Get output TTrees and add new branches to store matches
	TreeManager* muMgr = skroot_get_mgr(&muWriterLUN);
	TTree* muTree = muMgr->GetOTree();
	muTree->Branch("HwClockTicks", &HwClockTicks);
	muTree->Branch("NumRollovers", &NumRollovers);
	muTree->Branch("MatchedEvNums", &MatchedEvNums);
	muTree->Branch("MatchedEntryNums", &MatchedEntryNums);
	muTree->Branch("MatchedEntryHasAFT", &MatchedHasAFTs);
	muTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	muTree->Branch("MatchedParticleE", &MatchedParticleE);
	
	TreeManager* relicMgr = skroot_get_mgr(&relicWriterLUN);
	TTree* relicTree = relicMgr->GetOTree();
	relicTree->Branch("HwClockTicks", &HwClockTicks);
	relicTree->Branch("NumRollovers", &NumRollovers);
	relicTree->Branch("MatchedEvNums", &MatchedEvNums);
	relicTree->Branch("MatchedEntryNums", &MatchedEntryNums);
	relicTree->Branch("MatchedEntryHasAFT", &MatchedHasAFTs);
	relicTree->Branch("MatchedTimeDiff", &MatchedTimeDiff);
	relicTree->Branch("MatchedParticleE", &MatchedParticleE);
	
	return true;
}

bool ReconstructMatchedMuons::Execute(){
	
	Log(m_unique_name+" Relics to Write out: "+toString(m_data->writeOutRelics.size())+
	                  ", muons to write out: "+toString(m_data->muonsToRec.size()),v_debug,m_verbose);
	
	// write finished candidates to file
	if(m_data->writeOutRelics.size()){
		WriteEventsOut(m_data->writeOutRelics, relicWriterLUN, EventType::LowE);
	}
	
	if(m_data->muonsToRec.size()){
		WriteEventsOut(m_data->muonsToRec, muWriterLUN, EventType::Muon);
	}
	
}

bool ReconstructMatchedMuons::Finalise(){
	
	// write any remaining candidates to file
	if(m_data->writeOutRelics.size()){
		WriteEventsOut(m_data->writeOutRelics, relicWriterLUN, EventType::LowE);
	}
	
	if(m_data->muonsToRec.size()){
		WriteEventsOut(m_data->muonsToRec, muWriterLUN, EventType::Muon);
	}
	
	return true;
}


bool ReconstructMatchedMuons::WriteEventsOut(std::vector<ParticleCand>& eventsToWrite, int outLUN, EventType eventType){
	
	uint64_t currentEntry = rfmReader->GetEntryNumber();
	
	for(int i = 0; i < eventsToWrite.size(); i++){
		
		Log(m_unique_name+" Writing out next "+toString(eventType)+"; entry "+
		    toString(eventsToWrite[i].OutEntryNumber),v_debug,m_verbose);
		
		// first we need to (re-)load the muon/relic entry to save
		
		// N.B. we need to keep all hits (can't call delete_outside_hits), because:
		// For muons we need the entire mu+AFT window for the neutron cloud search
		// For relics we need the SHE+AFT window for the neutron tag search
		// In fact, we'll also merge the AFT hits with the primary readout here,
		// as dealing with pairs of TTree entries makes everything downstream more complicated.
		// In that case we need to re-load both entries, but for the AFT we're really only interested
		// in the hits; we'll keep the Header branch (event number, time etc) of the primary event.
		// In which case it's probably easier to read the AFT first, and just buffer the hits.
		static rawtqinfo_common rawtqinfo_aft;
		if(eventsToWrite[i].hasAFT && (eventsToWrite[i].InEntryNumber+1) != currentEntry){
			
			std::cout<<m_unique_name+" Rolling back input reader to grab AFT for "<<eventType
			         <<" entry "<<(eventsToWrite[i].InEntryNumber+1)<<std::endl;
			
			Log(m_unique_name+" prefetching AFT entry "+toString(eventsToWrite[i].InEntryNumber+1),
			    v_debug,m_verbose);
			
			// read the AFT entry from file
			get_ok = m_data->getTreeEntry(rfmReaderName, eventsToWrite[i].InEntryNumber+1);
			if(!get_ok){
				Log(m_unique_name+" Error reading AFT entry " 
				    +toString(eventsToWrite[i].InEntryNumber+1),v_error,m_verbose);
				return false;
			}
			
			// make a note of the AFT hits
			rawtqinfo_aft = rawtqinfo_;
			
		}
		
		// ok now grab the primary lowe/muon event
		
		// if we're going to be doing muon reconstruction we should reload
		// with noisy channels masked
		if(eventType==EventType::Muon){
			current_badch_masking = combad_.imaskbadopt;
			int newbadopt = 0;     // XXX 0 = mask all kinds of bad channels, NOT disable masking!
			skbadopt_(&newbadopt); // n.b. all this does is update the common block value
		}
		
		std::cout<<m_unique_name+" Rolling back reader to write out "<<eventType
		         <<" entry "<<eventsToWrite[i].InEntryNumber<<std::endl;
		
		m_data->getTreeEntry(rfmReaderName, eventsToWrite[i].InEntryNumber);
		
		// if this was a subtrigger we should also shift the time window accordingly
		if(eventsToWrite[i].SubTriggerNumber!=0){
			set_timing_gate_(&eventsToWrite[i].it0xsk);
			int neglun = -std::abs(rfmReaderLUN);
			skcread_(&neglun, &get_ok);
			// get_ok = 0 (physics entry), 1 (error), 2 (EOF), other (non-physics)
			if(get_ok!=0){
				Log(m_unique_name+" Error! skcread returned "+toString(get_ok)
				    +" when reloading subtrigger!",v_error,m_verbose);
				return false;
			}
		}
		
		if(eventType==EventType::Muon){
			// reset the bad channel masking for lowe events
			skbadopt_(&current_badch_masking);
		}
		
		if(eventType==EventType::Muon){
			// ok, do muon reconstruction
			// this will populate the reco_muons vector with a set of skroot_mu_common
			// objects, each representing the result of skroot_mu_ for each muon.
			// we need a vector as muboy can reconstruct multiple muons with unique entry points
			// and we calculate a separate dEdx array for each. (FIXME do we really need to do this?)
			get_ok = ReconstructNextMuon();
			if(!get_ok){
				Log(m_unique_name+" Error reconstructing muon event "
				    +toString(eventsToWrite[i].InEntryNumber),v_error,m_verbose);
				continue; // skip writing out this muon i guess...?
			}
		}
		
		// add those hits from the AFT to the primary event
		AddAftHits(rawtqinfo_aft);
		
		// set header and tq info
		// HEAD branch from assorted skhead_* common blocks,
		// TQREAL branch from rawtqinfo_ common block
		skroot_set_tree_(&outLUN);
		
		// update branch variables w/ info about matches
		MatchedEvNums = eventsToWrite[i].matchedParticleEvNum;
		MatchedEntryNums = eventsToWrite[i].matchedParticleEntryNum;
		MatchedHasAFTs = eventsToWrite[i].matchedParticleHasAFT;
		MatchedTimeDiff = eventsToWrite[i].matchedParticleTimeDiff;
		MatchedParticleE = eventsToWrite[i].matchedParticleBSEnergy;
		HwClockTicks = eventsToWrite[i].EventTicks;
		NumRollovers = eventsToWrite[i].NumRollovers;
		
		// for LowE events we need to set the LowE reconstruction info
		if(eventType==EventType::LowE){
			skroot_lowe_ = eventsToWrite[i].LowECommon;
			skroot_set_lowe_(&relicWriterLUN,
			                 skroot_lowe_.bsvertex,
			                 skroot_lowe_.bsresult,
			                 skroot_lowe_.bsdir,
			                 skroot_lowe_.bsgood,
			                 &skroot_lowe_.bsdirks,
			                 skroot_lowe_.bseffhit,
			                 &skroot_lowe_.bsenergy,
			                 &skroot_lowe_.bsn50,
			                 &skroot_lowe_.bscossun,
			                 skroot_lowe_.clvertex,
			                 skroot_lowe_.clresult,
			                 skroot_lowe_.cldir,
			                 &skroot_lowe_.clgoodness,
			                 &skroot_lowe_.cldirks,
			                 skroot_lowe_.cleffhit,
			                 &skroot_lowe_.clenergy,
			                 &skroot_lowe_.cln50,
			                 &skroot_lowe_.clcossun,
			                 &skroot_lowe_.latmnum,
			                 &skroot_lowe_.latmh,
			                 &skroot_lowe_.lmx24,
			                 &skroot_lowe_.ltimediff,
			                 &skroot_lowe_.lnsratio,
			                 skroot_lowe_.lsdir,
			                 &skroot_lowe_.spaevnum,
			                 &skroot_lowe_.spaloglike,
			                 &skroot_lowe_.sparesq,
			                 &skroot_lowe_.spadt,
			                 &skroot_lowe_.spadll,
			                 &skroot_lowe_.spadlt,
			                 &skroot_lowe_.spamuyn,
			                 &skroot_lowe_.spamugdn,
			                 skroot_lowe_.posmc,
			                 skroot_lowe_.dirmc,
			                 skroot_lowe_.pabsmc,
			                 skroot_lowe_.energymc,
			                 &skroot_lowe_.darkmc,
			                 &skroot_lowe_.islekeep,
			                 &skroot_lowe_.bspatlik,
			                 &skroot_lowe_.clpatlik,
			                 &skroot_lowe_.lwatert,
			                 &skroot_lowe_.lninfo,
			                 skroot_lowe_.linfo);
			
			// invoke TTree::Fill
			skroot_fill_tree_(&outLUN);
			
			/*
			// no longer needed as we've merged with the parent event
			// if the event had an associated AFT trigger, write that out as well
			if(eventsToWrite[i].hasAFT){
				Log(m_unique_name+" Writing out follow-up AFT entry",v_debug,m_verbose);
				m_data->getTreeEntry(rfmReaderName, eventsToWrite[i].InEntryNumber+1);
				skroot_set_tree_(&outLUN);
				skroot_fill_tree_(&outLUN);
			}
			*/
			
		}
		
		// for muon events, set the reconstructed muon info if available
		else if(eventType==EventType::Muon){
			
			std::vector<skroot_mu_common> mu_reco_info;
			get_ok = m_data->CStore.Get("reco_muons", mu_reco_info);
			
			// ok, so muboy may reconstruct multiple muons, and we should save all of them
			if(get_ok && mu_reco_info.size()){
				for(int j=0; j<mu_reco_info.size(); ++j){
					// note these only differ by the muon entry point
					// (for which we need to get the corresponding element from
					// skroot_mu_.muboy_entpos, using the index in skroot_mu_.mu_info[7])
					// and in the dE/dx arrays in mu_info[10:210] and skroot_mu_.muboy_dedx
					// to be honest we could save a lot of disk space by being smarter
					// in the way we read these in, but for now just make every entry independent.
					skroot_mu_ = mu_reco_info.at(j);
					
					skroot_set_mu_(&outLUN,
					               skroot_mu_.muentpoint,
					               skroot_mu_.mudir,
					               &skroot_mu_.mutimediff,
					               &skroot_mu_.mugoodness,
					               &skroot_mu_.muqismsk,
					               &skroot_mu_.muyn,
					               &skroot_mu_.mufast_flag,
					               &skroot_mu_.muboy_status,
					               &skroot_mu_.muboy_ntrack,
					               skroot_mu_.muboy_entpos,
					               skroot_mu_.muboy_dir,
					               &skroot_mu_.muboy_goodness,
					               &skroot_mu_.muboy_length,
					               skroot_mu_.muboy_dedx,
					               skroot_mu_.mubff_entpos,
					               skroot_mu_.mubff_dir,
					               &skroot_mu_.mubff_goodness,
					               &skroot_mu_.muninfo,
					               skroot_mu_.muinfo);
					
					// invoke TTree::Fill
					skroot_fill_tree_(&outLUN);
					
				}
			}
		}
		
	}
	
	eventsToWrite.clear();
	
	// reload the previous entry so that no issues are caused with other tools
	if(currentEntry != rfmReader->GetEntryNumber()){
		m_data->getTreeEntry(rfmReaderName, currentEntry);
	}
	
	return true;
}

bool ReconstructMatchedMuons::AddAftHits(const rawtqinfo_common& rawtqinfo_aft){
	
	// add hits from AFT tqinfo to currently loaded (prompt) event
	// skroot_fill_tree populates the TQREAL branch from the rawtqinfo_ common block;
	// in particular icabbf_raw, qbuf_raw, tbuf_raw for ID,
	// and icabaz_raw, qaskz_raw, taskz_raw for OD, so we need to copy over
	// the new entries for these arrays
	
	// we could think about skipping out-of-gate hits, since they're generally
	// not used for anything.
	
	// Han's ReadHits tool actually merges SHE+AFT and expects there to be
	// an overlap between the hit contents of the two,
	// so perhaps we should check for that and avoid duplicates...
	int k=0;
	int last_hit = rawtqinfo_.nqisk_raw-1;
	for(k=0; k<rawtqinfo_aft.nqisk_raw; ++k){
		if(rawtqinfo_aft.icabbf_raw[k] != rawtqinfo_.icabbf_raw[last_hit]) continue;
		if(rawtqinfo_aft.qbuf_raw[k]   != rawtqinfo_.qbuf_raw[last_hit])   continue;
		//if(rawtqinfo_aft.tbuf_raw[k]   != rawtqinfo_.tbuf_raw[last_hit])   continue;
		// can't compare times as they are relative to different trigger times
		break;
	}
	
	// i'm not aware of any particular sorting of hits required
	// so we'll just append the new hits to the end. Hope that works...!
	
	// ID hits
	int j = (k != rawtqinfo_aft.nqisk_raw) ?  k : rawtqinfo_.nqisk_raw;
	for(int i=0; i<rawtqinfo_aft.nqisk_raw; ++i){
		rawtqinfo_.tbuf_raw[j] = rawtqinfo_aft.tbuf_raw[i];
		rawtqinfo_.qbuf_raw[j] = rawtqinfo_aft.qbuf_raw[i];
		rawtqinfo_.icabbf_raw[j] = rawtqinfo_aft.icabbf_raw[i];
		//rawtqinfo_.itiskz_raw[j] = rawtqinfo_aft.itiskz_raw[i];
		//rawtqinfo_.iqiskz_raw[j] = rawtqinfo_aft.iqiskz_raw[i];
		// latter 2 not used by skroot_set_tree, so don't need to copy
		++j;
	}
	rawtqinfo_.nqisk_raw = j;
	
	// OD hits
	k=0;
	last_hit = rawtqinfo_.nhitaz_raw-1;
	for(k=0; k<rawtqinfo_aft.nhitaz_raw; ++k){
		if(rawtqinfo_aft.icabaz_raw[k] != rawtqinfo_.icabaz_raw[last_hit]) continue;
		if(rawtqinfo_aft.qaskz_raw[k]   != rawtqinfo_.qaskz_raw[last_hit]) continue;
		//if(rawtqinfo_aft.taskz_raw[k]   != rawtqinfo_.taskz_raw[last_hit]) continue;
		// can't compare times as they are relative to different trigger times
		break;
	}
	
	j = (k != rawtqinfo_.nhitaz_raw) ? k : rawtqinfo_.nhitaz_raw;
	for(int i=0; i<rawtqinfo_aft.nhitaz_raw; ++i){
		rawtqinfo_.taskz_raw[j] = rawtqinfo_aft.taskz_raw[i];
		rawtqinfo_.qaskz_raw[j] = rawtqinfo_aft.qaskz_raw[i];
		rawtqinfo_.icabaz_raw[j] = rawtqinfo_aft.icabaz_raw[i];
		//rawtqinfo_.itaskz_raw[j] = rawtqinfo_aft.itaskz_raw[i];
		//rawtqinfo_.iqaskz_raw[j] = rawtqinfo_aft.iqaskz_raw[i];
		++j;
	}
	rawtqinfo_.nhitaz_raw = j;
	
	return true;
}


bool ReconstructMatchedMuons::ReconstructNextMuon(){
	
	// we'll return a vector of skroot_mu_common structs.
	// why a vector? because muboy might reconstruct multiple muons.
	reco_muons.clear();
	
	// for relic spallation checks we only bother with BFF for the muon
	// if the corresponding relic bsenergy is > 12... not sure the logic of that,
	// but we may wish to do a similar thing.
	bool tryBFF=false;
	m_data->vars.Get("tryBFF",tryBFF);
	
	// =================== //
	// Muon Reconstruction //
	// =================== //
	
	// store charge ranges before fix_maxqisk
	skroot_mu_.muinfo[0] = skq_.qismsk;
	skroot_mu_.muinfo[2] = skq_.qimxsk;
	
	// supposedly this undoes an upstream charge saturation correction
	// which was required for SKI-III but is no longer applicable for SKIV+
	fix_maxqisk_();
	
	// save updated charge metrics
	skroot_mu_.muqismsk = skq_.qismsk;
	skroot_mu_.muinfo[3] = skq_.qimxsk;  // should we not update the value in skq_.?
	if(skroot_mu_.muninfo < 4) skroot_mu_.muninfo = 4;
	
	int muyn_org, muynf;
	
	//Muon reconstruction developed by Tomoeda and Yamaguchi
	mfmuselect_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &skroot_mu_.mugoodness, &muyn_org);
	
	//muyn == 1 - good fit
	//muyn == 0 - bad fit
	
	if(muyn_org > 0){
		skroot_mu_.muyn = 1;
	} else if(muyn_org < 0) {
		skroot_mu_.muyn = 0;
	} else {
		Log("Muyn_org returning as == 0. Not supported yet", v_error, m_verbose);
		return false;
	}
	
	//Apply fast fit if mfmuselect has returned a bad fit
	if(skroot_mu_.muyn == 0){
		mffastfast_(&skroot_mu_.muentpoint, &skroot_mu_.mudir, &muynf);
		skroot_mu_.mufast_flag = 1;
	}else{
		skroot_mu_.mufast_flag = 0;
	}
	
	skroot_mu_.muyn = muyn_org;
	if(skroot_mu_.muyn == 0){
		skroot_mu_.muyn = muynf;
		// XXX n.b. mufit_sk4_loose overwrites this with 'calflag', flagging various types of trigger..?
	}
	
	//Apply muboy
	int n_left;
	float  muentry[4], muboy_otherentry[36];
	// $ATMPD_ROOT/src/recon/fit/muboy.F
	muboy_zbs_(&skhead_.nevsk,
	           &skroot_mu_.muboy_status,    // stopping, through-going, corner-clipper, etc. 0=fit failed.
	           &muentry,                    // [0-2]: pos of PMT closest to entry point?, [3]: entry time
	           &skroot_mu_.muboy_dir,       // primary direction at entry point, unit normalised
	           &skroot_mu_.muboy_length,    // track length [cm]
	           &skroot_mu_.muboy_goodness,  // 0-1, higher is better*
	           &skroot_mu_.muboy_ntrack,    // num tracks ("can be 1 if multiple muons" 🤦)
	           &muboy_otherentry,           // additional entry points for tracks 2-5
	           &n_left);                    // num hit PMTs left after cluster finding
	
	// get muon track entry position(s)
	for(int track = 0; track < skroot_mu_.muboy_ntrack; track++){
		if(track == 0){
			skroot_mu_.muboy_entpos[track][0] = muentry[0];
			skroot_mu_.muboy_entpos[track][1] = muentry[1];
			skroot_mu_.muboy_entpos[track][2] = muentry[2];
			skroot_mu_.muboy_entpos[track][3] = muentry[3];
		}else{
			skroot_mu_.muboy_entpos[track][0] = muboy_otherentry[4*track - 7];
			skroot_mu_.muboy_entpos[track][1] = muboy_otherentry[4*track - 6];
			skroot_mu_.muboy_entpos[track][2] = muboy_otherentry[4*track - 5];
			skroot_mu_.muboy_entpos[track][3] = muboy_otherentry[4*track - 4];
		}
	}
	
	if(m_verbose > v_debug+2){
		std::cout << "muboy result:\n"
		          <<"\tgoodness: "<< skroot_mu_.muboy_goodness << "\n"
		          <<"\tntracks: " << skroot_mu_.muboy_ntrack   << "\n"
		          << "\tclass: "  << skroot_mu_.muboy_status   << "\n"
		          << "\tdir: ("   << skroot_mu_.muboy_dir[0]   << ", "
		                          << skroot_mu_.muboy_dir[1]   << ", "
		                          << skroot_mu_.muboy_dir[2]   << ")\n"
		          << "\tlength: " << skroot_mu_.muboy_length   << "\n"
		          << "\tfirst track entry point: ("
		                          << skroot_mu_.muboy_entpos[0][0] << ", "
		                          << skroot_mu_.muboy_entpos[0][1] << ", "
		                          << skroot_mu_.muboy_entpos[0][2] << ")\n"
		          <<"\ttime: "    << skroot_mu_.muboy_entpos[0][3] << "\n"
		          << std::endl;
	}
	
	// makededx needs primary entry position and direction.
	// copy muboy values to a set of variables which may overridden by BFF
	// (we already have muentry for position)
	float mudir[3];
	for(int i=0; i<4; ++i){
		if(i<3) mudir[i] = skroot_mu_.muboy_dir[i];
		muentry[i] = skroot_mu_.muboy_entpos[i][0];
	}
	
	// according to Scott's (ambiguously worded) lowe school slide:
	// for single through-going muons, stopping muons (status 1,2 respectively),
	// or sometimes for large showers, a goodness of 0.4-0.6+ is a good fit?
	// goodness values down to 0.15 are sometimes ok
	// goodness <0.1 are always bad.
	// (think that's the right interpretation of his slide?)
	
	// if muboy failed, try BFF... but as this can take up to 30 minutes per muon(!),
	// we have two flags: one Tool config which vetos the use of BFF on anything,
	// and one in `m_data->vars` which can veto based on the details of the event.
	if(!noBFF && tryBFF && skroot_mu_.muboy_status == 1 && skroot_mu_.muboy_goodness < 0.4){
		
		Log(m_unique_name+": Muboy failed, trying BFF",v_error,m_verbose);
		
		float bffpos[3];
		float hpos[3];
		newmufit_(&bffpos, &hpos, &skroot_mu_.mubff_goodness);
		
		Log(m_unique_name+": Finished BFF",v_error,m_verbose);
		
		// copy out result
		float modd = sqrt( pow((hpos[0]-bffpos[0]),2) + pow((hpos[1]-bffpos[1]),2) + pow((hpos[2]-bffpos[2]),2) );
		for(int j = 0; j < 3; j++){
			skroot_mu_.mubff_entpos[j] = bffpos[j];
			skroot_mu_.mubff_dir[j] = (hpos[j] - bffpos[j])/modd;
		}
		
		// if BFF succeeded, update the primary muon entry position and direction for makededx
		if(skroot_mu_.mubff_goodness > 0.3){
			for(int j = 0; j < 3; j++){
				muentry[j] = skroot_mu_.mubff_entpos[j];
				mudir[j] = skroot_mu_.mubff_dir[j];
			}
			
			// * muinfo[6]: flag for whether dedx was calculated using BFF or muboy.
			skroot_mu_.muinfo[6] = 1;
			
		} else {
			skroot_mu_.muinfo[6] = 0;
		}
		
	} else {
		// otherwise initialise BFF
		skroot_mu_.mubff_goodness = 0;
		for(int i=0; i<3; ++i){
			skroot_mu_.mubff_entpos[i] = 0;
			skroot_mu_.mubff_dir[i] = 0;
		}
		skroot_mu_.muinfo[6] = 0;
	}
	
	// mufit_sk4_loose saves all muboy tracks as separate events,
	// recalculating dedx for each with the corresponding muboy entry point.
	// XXX muboy tracks >1 are saved even when the primary track is overwritten by BFF
	// (although perhaps the muboy_status==1 BFF check implies muboy_ntrack==1)
	// ughhhhhhhhhhh we're gonna need another loop to save these
	for(int mutrack=0; mutrack<skroot_mu_.muboy_ntrack; ++mutrack){
		
		skroot_mu_.muinfo[7] = mutrack;
		
		// XXX note muboy only provides one direction, even for multiple muons.
		if(mutrack>0){
			muentry[0] = skroot_mu_.muboy_entpos[mutrack][0];
			muentry[1] = skroot_mu_.muboy_entpos[mutrack][1];
			muentry[2] = skroot_mu_.muboy_entpos[mutrack][2];
			muentry[3] = skroot_mu_.muboy_entpos[mutrack][3];
		}
		
		// calculate rate of energy loss along track
		/*
		// from $ATMPD_ROOT/src/recon/stmu/makededx.F
		// called by lomu_gd.F, which makes lomu_gd files
		makededx_(&muentry,
		          &mudir,
		          &skchnl_.ihcab,
		          &skq_.qisk,
		          &skt_.tisk,
		          &geopmt_.xyzpm,
		          &skq_.nqisk,
		          &skroot_mu_.muboy_dedx);  // populates this
		// n.b. there is also an updated version: $ATMPD_ROOT/src/recon/stmu/makededx2.F
		// that makes water transparency and coverage corrections...
		*/
		
		// but apparently this (either?) dedx is somehow 'not good'...?
		// the previous SRN analyses recalculate it with two other methods: kirk's dedx and scott's dedx.
		// FIXME understand what the issue is here and if, or why, we need to call all of them...
		
		// Kirk's method.
		// the relic_sk4_ana repo slightly modified kirk's dedx, accepting water transparency
		// and applying it in an expl function. In kirk's original version the water transparency
		// is looked up internally with lfwatert, but is never used....
		/*
		float watert;
		get_ok = m_data->vars.Get("watert",watert);
		if(!get_ok){
			// the TreeReader should maintain a suitable value in m_data->vars...
			// is the file MC, but no reference run was provided?
			Log(m_unique_name+" Error! relic_sk4_ana's `makededx` requires water transparency, "
				"but none found in m_data->vars!",v_error,m_verbose);
			return false;
		}
		*/
		
		// recalculate dedx with kirk's method
		// we'll put kirk's dedx array (200 elements) in muinfo[10]:muinfo[210]
		// previous relic sk4 code saved it in muinfo from element 0, but that overwrites:
		// (from $SKOFL_ROOT/inc/lowe/skroot_loweC.h)
		// * muinfo[0]: qismsk before fix_maxqisk
		// * muinfo[1]: subtrigger number
		// * muinfo[2]: original qimxsk before fix_maxqisk
		// * mufino[3]: qimxsk after fix_maxqisk (redundant with tq_.qimxsk?)
		// * muinfo[4]: parent muon event number (parent muon event for a muon event?)
		// * muinfo[5]: subtrigger number in AFT (how does this differ from rmuinfo[1]?)
		// we also make use of a couple more:
		// * muinfo[6]: whether BFF was applied (1) or not (0)
		// * muinfo[7]: muboy track number (for distinguishing later muboy tracks)
		// i can't find anything that seems to use the remaining elements, so just start from 10.
		float (*muinfo)[200] = (float(*)[200])(skroot_mu_.muinfo+10);
		// from $RELIC_WORK_DIR/lomufit/{lowfit/mufit}/src/makededx.F
		makededx_(&muentry,
		          &mudir,
		          &skchnl_.ihcab,
		          &skq_.qisk,
		          &skt_.tisk,
		          &geopmt_.xyzpm,
		          &skq_.nqisk,
		          &skhead_.nrunsk,
		//        &watert,*
		          muinfo);       // populates this
		// * uncomment this, and the appropriate signature in fortran_routines.h to use relic_sk4_ana ver
		// TODO maybe we could put them in namespaces to avoid conflicting names
		
		// skroot_get_mu_ initialises to 0 all muinfo elements > muninfo
		skroot_mu_.muninfo = 210;
		
		// Scott's method.
		// FIXME seems like we're just overwriting the original results with this.
		// if we're just going to overwrite the result, why even call the original makededx?
		// still, this is what lomufit_gd does for the "official" lomugd files,
		// so presumably there's some reason. Maybe the original value is a prerequisite???
		// from $SKOFL_ROOT/lowe/sklowe/makededx_intg.cc
		makededx_intg_(&muentry,
		               &mudir,
		               &skroot_mu_.muboy_length,
		               &skchnl_.ihcab,
		               &skq_.qisk,
		               &skt_.tisk,
		               &geopmt_.xyzpm,
		               &sktqz_.nqiskz,
		               &skhead_.nrunsk,
		               &skroot_mu_.muboy_dedx,    // populates this
		               &sktqz_.ihtiflz,
		               &skhead_.nevsk);
		
		/* can't do this if we have multiple muons reconstructed per Execute...
		// pass reconstructed muon info to output Tree branch variables
		skroot_set_mu_(&lun,
		               skroot_mu_.muentpoint,
		               skroot_mu_.mudir,
		               &skroot_mu_.mutimediff,
		               &skroot_mu_.mugoodness,
		               &skroot_mu_.muqismsk,
		               &skroot_mu_.muyn,
		               &skroot_mu_.mufast_flag,
		               &skroot_mu_.muboy_status,
		               &skroot_mu_.muboy_ntrack,
		               skroot_mu_.muboy_entpos,
		               skroot_mu_.muboy_dir,
		               &skroot_mu_.muboy_goodness,
		               &skroot_mu_.muboy_length,
		               skroot_mu_.muboy_dedx,
		               skroot_mu_.mubff_entpos,
		               skroot_mu_.mubff_dir,
		               &skroot_mu_.mubff_goodness,
		               &skroot_mu_.muninfo,
		               skroot_mu_.muinfo);
		*/
		reco_muons.push_back(skroot_mu_);
		
	}
	// ======================= //
	// End Muon Reconstruction //
	// ======================= //
	
	return true;
}


