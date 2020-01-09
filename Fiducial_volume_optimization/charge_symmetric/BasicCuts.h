// --- Define basic cuts here ---
double x=3300; //Default BiasFlashtime variable
TCut cRandom("EventCategory == 1");
TCut cNotRandom("EventCategory != 1");
TCut cNotEmpty("Empty == 0");
TCut cGoodFlash("BiasFlashTime < 10800"); //10800 for T5Z2 // 3300 for T2Z1// 3000 for Ylite

double flashtime()
{
return x;
}

int pstd=0; // default choose T2Z1 for PStd variable
int bt=0;  // default choose T2Z1 for bias flashtime variable

TCut cBaseTemp("BaseTemp>0");

TCut cVoltageBiasS1("QIS1bias>1.95 && QIS1bias<2.05 && QOS1bias>1.95 && QOS1bias<2.05");
TCut cVoltageBiasS2("QIS2bias<-1.95 && QIS2bias>-2.05 && QOS2bias<-1.95 && QOS2bias>-2.05");
TCut cVoltageBias=cVoltageBiasS1+cVoltageBiasS2;

TCut cGoodStartTime("PTOFdelay > -195e-6 && PTOFdelay < 35e-6");

TCut ch3Cuts=cNotRandom+cGoodFlash;
TCut cBasicCuts=cNotEmpty+cGoodFlash+cBaseTemp+cVoltageBias+cGoodStartTime;
//TCut cBasicCuts2=cRandom+cNotEmpty+cGoodFlash+cBaseTemp+cPstd+cVoltageBias;
//-----------------------------------------For CDMSlite-------------------------------------
/*
//TCut cHVvolts("HVvolts < -68 && HVvolts > -72");
TCut cHVvolts("HVvolts < -23 && HVvolts > -27");
TCut cLeakCurr("HVnamps<=3.0");
TCut cZeroCurr("HVnamps>0.0");
TCut cBaseTemplite("BaseTemp>0");
TCut cGoodStartTime("PTOFdelay > -195e-6 && PTOFdelay < 35e-6"); // zip14
//TCut cGoodStartTime("PTOFdelay > -200e-6 && PTOFdelay < 35e-6"); //zip4
TCut cBasicCutslite=cNotRandom+cNotEmpty+cGoodFlash+cHVvolts+cZeroCurr+cBaseTemplite+cLeakCurr+cGoodStartTime;

//TCut cLEChisq("PTNFchisq < (0.015*ptNF*ptNF+4600) && PTNFchisq >3600");
TCut cLEChisq("PTNFchisq < (0.015*ptNF*ptNF+4600)");
TCut cDeltaChisqGlitch("(PTOFchisq-PTglitch1OFchisq)<((-2.3*(ptNF*ptNF))+25)");
TCut cDeltaChiSqLFN("(PTOFchisq-PTlfnoise1OFchisq)<((-4.54*(ptNF*ptNF))+25)");
TCut cBadSeries("SeriesNumber!=11509181456 && SeriesNumber!=11506081315 && SeriesNumber!=11506081424 && SeriesNumber!=11506081534 && SeriesNumber!=11506300206");
TCut cQualityCutslite=cLEChisq+cDeltaChisqGlitch+cDeltaChiSqLFN+cBadSeries;*/
