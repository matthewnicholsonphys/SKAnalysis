#include "Factory.h"

Tool* Factory(std::string tool){
Tool* ret=0;

// if (tool=="Type") tool=new Type;
if (tool=="DummyTool") ret=new DummyTool;
if (tool=="ApplyTMVA") ret=new ApplyTMVA;
if (tool=="ExtractFeatures") ret=new ExtractFeatures;
if (tool=="ReadHits") ret=new ReadHits;
if (tool=="ReadMCInfo") ret=new ReadMCInfo;
if (tool=="SearchCandidates") ret=new SearchCandidates;
if (tool=="SetPromptVertex") ret=new SetPromptVertex;
if (tool=="SKRead") ret=new SKRead;
if (tool=="SubtractToF") ret=new SubtractToF;
if (tool=="NTupleMatcher") ret = new NTupleMatcher;
if (tool=="AddNoise") ret = new AddNoise;
if (tool=="WriteOutput") ret=new WriteOutput;
if (tool=="TestTool") ret=new TestTool;
if (tool=="TruthNeutronCaptures") ret=new TruthNeutronCaptures;
if (tool=="TruthNeutronCaptures_v2") ret=new TruthNeutronCaptures_v2;
if (tool=="TruthNeutronCaptures_v3") ret=new TruthNeutronCaptures_v3;
if (tool=="LoadFileList") ret=new LoadFileList;
if (tool=="RootReadTest") ret=new RootReadTest;
if (tool=="MakeNCaptTree") ret=new MakeNCaptTree;
if (tool=="GracefulStop") ret=new GracefulStop;
if (tool=="LoadBetaSpectraFluka") ret=new LoadBetaSpectraFluka;
if (tool=="PlotMuonDtDlt") ret=new PlotMuonDtDlt;
if (tool=="FitLi9Lifetime") ret=new FitLi9Lifetime;
if (tool=="FitPurewaterLi9NcaptureDt") ret=new FitPurewaterLi9NcaptureDt;
if (tool=="FitSpallationDt") ret=new FitSpallationDt;
if (tool=="PurewaterSpallAbundanceCuts") ret=new PurewaterSpallAbundanceCuts;
if (tool=="TreeReader") ret=new TreeReader;
if (tool=="lfallfit_simple") ret=new lfallfit_simple;
if (tool=="TreeReaderDemo") ret=new TreeReaderDemo;
if (tool=="lfallfit") ret=new lfallfit;
if (tool=="evDisp") ret=new evDisp;
if (tool=="CompareRootFiles") ret=new CompareRootFiles;
if (tool=="vectgen") ret=new vectgen;
if (tool=="SK2p2MeV_ntag") ret=new SK2p2MeV_ntag;
if (tool=="ntag_BDT") ret=new ntag_BDT;
if (tool=="SimplifyTree") ret=new SimplifyTree;
if (tool=="NCaptInfo_NTag") ret=new NCaptInfo_NTag;
if (tool=="NCaptInfo_BDT") ret=new NCaptInfo_BDT;
if (tool=="TrueNCaptures") ret=new TrueNCaptures;
if (tool=="ReadMCParticles") ret=new ReadMCParticles;
if (tool=="PlotNeutronCaptures") ret=new PlotNeutronCaptures;
if (tool=="PrintEvent") ret=new PrintEvent;
if (tool=="PlotHitTimes") ret=new PlotHitTimes;
if (tool=="VertexFitter") ret=new VertexFitter;
if (tool=="CombinedFitter") ret=new CombinedFitter;
if (tool=="MuonSearch") ret=new MuonSearch;
if (tool=="IDHitsCut") ret=new IDHitsCut;
if (tool=="ODCut") ret=new ODCut;
if (tool=="RunwiseEnergyCut") ret=new RunwiseEnergyCut;
if (tool=="WallCut") ret=new WallCut;
if (tool=="RelicCandidates") ret=new RelicCandidates;
if (tool=="RelicMuonMatching") ret=new RelicMuonMatching;
if (tool=="WriteSpallCand") ret=new WriteSpallCand;
if (tool=="ReconstructMatchedMuons") ret=new ReconstructMatchedMuons;

if (tool=="SkipTriggers") ret=new SkipTriggers;
if (tool=="BuildHist") ret=new BuildHist;
if (tool=="CutRecorder") ret=new CutRecorder;
if (tool=="SkipEventFlags") ret=new SkipEventFlags;
if (tool=="AddTree") ret=new AddTree;
if (tool=="WriteSkEvent") ret=new WriteSkEvent;
if (tool=="RelicMuonPlots") ret=new RelicMuonPlots;
return ret;
}

