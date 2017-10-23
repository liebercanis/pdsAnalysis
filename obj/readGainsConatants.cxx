// obsolete
Int_t pmtAna::readGainConstants(TString fileName)
{
  ifstream in;
  in.open(fileName);
  Int_t ngains=0;
  if(!in.is_open() ) {
    printf(" cannot open file %s \n",fileName.Data());
    return ngains;
  }
  Double_t rgain[NPMT];
  printf(" readGainConstants from file %s \n",fileName.Data());
  string line,label,type,sgain;
  while (in.good()) {
    in >> label >> type >> sgain;
    if(in.eof()) break;
    // look for comment or blank line
    if( label.find("%") != std::string::npos || label.size()<4 ) {
      getline(in,line); // throw away line
      continue;
    }
    if( label.size()<2) continue;
    //cout << label << "  " << type << "  " << sgain << endl;
    int b = atoi(&label[1]);
    int c = atoi(&label[3]);
    int ipmt = toPmtNumber(b,c);
    rgain[ipmt] = atof(sgain.c_str());
    ++ngains;
  }
  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" %i  %f \n",ipmt,rgain[ipmt]); 
  // normalize 
  for(int ipmt=0; ipmt<NPMT; ++ipmt) { gain[ipmt] = rgain[ipmt]/rgain[0];}
  printf(" normalized \n");
  for(int ipmt=0; ipmt<NPMT; ++ipmt) printf(" %i  %f ; ",ipmt,gain[ipmt]); 
  printf(" \n");

  Double_t gsumAve[NPMT];
  gsumAve[0]= 1.000 ;
  gsumAve[1]= 0.962 ;
  gsumAve[2]= 1.089 ;
  gsumAve[3]= 1.106 ;
  gsumAve[4]= 0.682 ;
  gsumAve[5]= 0.900 ;
  gsumAve[6]= 0.973 ;
  gsumAve[7]= 1.021 ;
  gsumAve[8]= 0.583 ;
  gsumAve[9]= 0.623 ;
  gsumAve[10]= 0.636 ;
  gsumAve[11]= 0.743 ;
  gsumAve[12]= 0.609 ;
  gsumAve[13]= 0.831 ;
  gsumAve[14]= 0.773 ;
  gsumAve[15]= 0.954 ;
  gsumAve[16]= 0.937 ;
  gsumAve[17]= 0.872 ;
  gsumAve[18]= 1.225 ;
  gsumAve[19]= 1.206 ;
  gsumAve[20]= 0.858 ;

 
  printf(" my gains \n");

  for(int ipmt=0; ipmt<NPMT; ++ipmt) {
    gain[ipmt]=1;
    printf(" %i  %f  ; ",ipmt,gain[ipmt]); 
  }

  printf(" \n");

  in.close();
  return ngains;
}


