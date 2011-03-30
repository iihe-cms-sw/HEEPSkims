{

gROOT->ProcessLine(".L DataMCTreeAnalysis.C+")   ;
a = new DataMCTreeAnalysis();
a->Loop()  ;

}
