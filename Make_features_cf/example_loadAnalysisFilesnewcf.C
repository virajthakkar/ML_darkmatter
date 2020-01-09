//////////////////////////////////////////////////////////////////
//
// This is a set of utility functions that will read specified data
// and then friend all the relevant trees
//
// To use it:
//
// 1. Specify the path to the rq/rrq/cut files and pass it into
// chainDataAll as the dataDir argument.
//
// Original author of this script:  Scott Fallows (some mods by LLH)
//
//////////////////////////////////////////////////////////////////


//chains and friends all files, including cuts
TChain* chainDataAllSpecial(Int_t z, string dataDir, string week, string type, Int_t z2) 
{

   if(dataDir == "")
   {
     cout <<"ERROR! Environmental variable CDMS_ANALYSISDATA is not set!" << endl;
   }

   string base = dataDir;
   cout <<"Reading data from viraj's directory: \n" << base << endl;

   TChain* eventChain = new TChain( "eventTree" );
   TChain* mergeChain = new TChain( Form("zip%d",z) );
   TChain* calibChain = new TChain( Form("calibzip%d",z) );
   TChain* calibEvent = new TChain( "calibevent" );

   
/*   eventChain->Add( Form("%s/merge_Prodv5-6-3_%s*_%s_z%dlite*.root/rqDir/eventTree",base.c_str(), week.c_str(), type.c_str(), z2) );
   mergeChain->Add( Form("%s/merge_Prodv5-6-3_%s*_%s_z%dlite*.root/rqDir/zip%d",base.c_str(), week.c_str(), type.c_str(), z2, z) );
   calibChain->Add( Form("%s/calib_Prodv5-6-3_%s*_%s_z%dlite*.root/rrqDir/calibzip%d",base.c_str(), week.c_str(), type.c_str(), z2,z) );
   calibEvent->Add( Form("%s/calib_Prodv5-6-3_%s*_%s_z%dlite*.root/rrqDir/calibevent",base.c_str(), week.c_str(), type.c_str(), z2) );
   // typical name: merge_Prodv5-6-3_1410a_cf_z14lite.root   
*/
   eventChain->Add( Form("%s/merge_Prodv5-3-5_z%dlite_%s*_%s.root/rqDir/eventTree",base.c_str(),z2,week.c_str(),type.c_str()) );
   mergeChain->Add( Form("%s/merge_Prodv5-3-5_z%dlite_%s*_%s.root/rqDir/zip%d",base.c_str(), z2,week.c_str(),type.c_str(), z) );
   calibChain->Add( Form("%s/calib_Prodv5-3-5_z%dlite_%s*_%s.root/rrqDir/calibzip%d",base.c_str(), z2, week.c_str(),type.c_str(), z) );
   calibEvent->Add( Form("%s/calib_Prodv5-3-5_z%dlite_%s*_%s.root/rrqDir/calibevent",base.c_str(), z2, week.c_str(),type.c_str()) );
// typical name : merge_Prodv5-3-5_z14lite_1401_ba.root
   cout<<"eventchain"<<endl;
cout<<Form("%s/merge_Prodv5-3-5_z%dlite_%s*_%s.root/rqDir/eventTree",base.c_str(),z2,week.c_str(),type.c_str());

cout<<"mergeChain"<<endl;
cout<<Form("%s/merge_Prodv5-3-5_z%dlite_%s*_%s.root/rqDir/zip%d",base.c_str(), z2,week.c_str(),type.c_str(), z) ;
//After chaining files together, friend trees
   
   if(eventChain->GetEntries() == mergeChain->GetEntries())
   {
     mergeChain->AddFriend(eventChain);
   }
   

   if(calibEvent->GetEntries() == mergeChain->GetEntries())
   {
     mergeChain->AddFriend(calibEvent);
   }
   

   if(calibChain->GetEntries() == mergeChain->GetEntries())
   {
     mergeChain->AddFriend(calibChain);
   }
   

   return mergeChain;
}

/*
//chains and friends all files, including cuts
TChain* chainEventOnly(string type="bg", string dataDir, Bool_t withcuts=1) 
{

  cout <<"Hello!" << endl;

   if(dataDir == "")
   {
     cout <<"ERROR! Environmental variable CDMS_ANALYSISDATA is not set!" << endl;
   }

   string base = dataDir + "/merged";
   cout <<"Reading data from directory: \n" << base << endl;

   TChain* eventChain = new TChain( "eventTree" );
   eventChain->Add( Form("%s/all/%s/merge*.root/rqDir/eventTree",base.c_str(),type.c_str()) );

   //retrieve, chain and friend cuts
   if( withcuts ) 
   {
      vector<string> cutlist;
      vector<Int_t> isgeneral;
      string acut;
      cutlist.clear();
      isgeneral.clear();
      
      //read the cut file
      string cutfilename = type + "cuts.txt";
      ifstream infile (cutfilename.c_str(), ios_base::in);
      if(infile.fail()) 
	{ 
	  cout <<"Error reading cut file !" << endl; 
	  exit();
	}
      
      //store cuts and their flags in a vector
      while (getline(infile, acut, '\n'))
      {
         cutlist.push_back( acut.substr(0,acut.size()-2) );
         isgeneral.push_back( atoi( (acut.substr(acut.size()-1,1)).c_str() ) );
      }

      //loop over cuts and for each, chain, then friend
      for(UInt_t cn=0; cn < cutlist.size(); cn++)
      {
         Char_t* cutName = (Char_t*)cutlist[cn].c_str();
	 
	 //skip if its a detector specific cut
	 if(isgeneral[cn] == 0)  continue;

	 int z = 0; //dummy zip value
	 TChain* cutChain = getCutChain(z, base.c_str(), type.c_str(), cutName, isgeneral[cn]); 

	 //chain if lengths match
	 if(cutChain->GetEntries() == eventChain->GetEntries())
	 {
	   eventChain->AddFriend(cutChain);
	 }
	 else
	 {
	   cout <<"ERROR! cutTree and eventTree lengths do not match!" 
		<<"\nCut tree length = " <<cutChain->GetEntries() << endl;
	   exit(1);
	 }
      }
   }

   return eventChain;

}
*/

