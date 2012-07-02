#include <boost/algorithm/string.hpp>
#include "msgfdb2Pin.h"
#include <ssl.h>

msgfdb2Pin::msgfdb2Pin() {
  reader = new msgfdbReader(); //Change this to a variable that gets inialised in run?
}

msgfdb2Pin::~msgfdb2Pin() {
  //deleting temporary folder(s)
  //NOTE destroy reader??
  if(reader)
    delete reader;
}

string msgfdb2Pin::extendedGreeter() {
  ostringstream oss;
  char* host = getenv("HOSTNAME");
  oss << greeter();
  oss << "Issued command:" << endl << call << endl;
  oss.seekp(-1, ios_base::cur);
  if(host) oss << "on " << host << endl;
  return oss.str();
}

string msgfdb2Pin::greeter() {
  ostringstream oss;
  oss << "\nmsgfdb2pin version " << VERSION << ", ";
  oss << "Build Date " << __DATE__ << " " << __TIME__ << endl;
  oss << "Copyright (c) 2010 Lukas Käll. All rights reserved." << endl;
  oss << "Written by Lukas Käll (lukask@cbr.su.se) in the" << endl;
  oss << "Department of Biochemistry and Biophysics at the Stockholm University."
      << endl;
  return oss.str();
}

bool msgfdb2Pin::parseOpt(int argc, char **argv) {
  ostringstream callStream;
  callStream << argv[0];
  for (int i = 1; i < argc; i++) {
    callStream << " " << argv[i];
  }
  callStream << endl;
  call = callStream.str();
  ostringstream intro, endnote;
  intro << greeter() << endl << "Usage:" << endl;
  intro << "   msgfdb2pin [options] -o output.xml target decoy " << endl << endl;
  intro << "Where output.xml is where the output will be written (ensure to have read and " << endl;
  intro << "write access on the file).target is the msgfdb file with the targets, and decoy is" << endl;
  intro << "the msgfdb file with the decoys. Small data sets may be merged by replace the files with" << endl;
  intro << "meta files. Meta files are text files containing the paths of files, one" << endl;
  intro << "path per line. For successful result, the different runs should be generated" << endl;
  intro << "under similar condition." << endl;

  // init
  CommandLineParser cmd(intro.str());

  cmd.defineOption("o",
      "outputXML",
      "save output in an XML file",
      "filename");
  cmd.defineOption("m",
      "matches",
      "Maximal number of matches to take in consideration per spectrum when using msgfdb-files",
      "number");
  cmd.defineOption("v",
      "verbose",
      "Set verbosity of output: 0=no processing info, 5=all, default is 2",
      "level");
  cmd.defineOption("u",
      "unitnorm",
      "Use unit normalization [0-1] instead of standard deviation normalization",		//NOTE not used?
      "",
      TRUE_IF_SET);
  cmd.defineOption("a",
      "aa-freq",
      "Calculate amino acid frequency features",
      "",
      TRUE_IF_SET);
  cmd.defineOption("b",
      "PTM",
      "Calculate feature for number of post-translational modifications",
      "",
      TRUE_IF_SET);
  cmd.defineOption("Q",
      "quadratic",
      "Calculate quadratic feature terms",							//?
      "",
      TRUE_IF_SET);
  cmd.defineOption("y",
      "notryptic",
      "Turn off calculation of tryptic/chymo-tryptic features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("c",
      "chymo",
      "Replace tryptic features with chymo-tryptic features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("t",
      "thermolysin",
      "Replace tryptic features with thermlysinic features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("k",
      "proteinasek",
      "Replace tryptic features with proteinase-k features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("n",
      "pepsin",
      "Replace tryptic features with peptic features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("e",
      "elastase",
      "Replace tryptic features with elastase features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("l",
      "lys-n",
      "Replace tryptic features with Lys-N features.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("N",
      "PNGaseF",
      "Calculate feature based on N-linked glycosylation pattern resulting from a PNGaseF treatment. (N[*].[ST])",
      "",
      TRUE_IF_SET);
  cmd.defineOption("2",
      "ms2-file",
      "File containing spectra and retention time. The file could be in mzXML, MS2 or compressed MS2 file.",
      "filename");
  cmd.defineOption("M",
      "isotope",
      "Mass difference calculated to closest isotope mass rather than to the average mass.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("p",
		   "psm-annotation",
      "An anotation scheme used to convert the psms from the search. An example if Q# was used to describe pyro-glu formation (UNIMOD:28), and S* and T* was used to describe phosphorylation (UNIMOD:21), we would use the option -p *:21:#:28",
      "Scheme");													//?
  cmd.defineOption("Z",
      "combined",
      "The input will be a target/decoy combined file.",
      "",
      TRUE_IF_SET);
  cmd.defineOption("P",
      "pattern",
      "Pattern used to identify the decoy PSMs",
      "",
      "pattern");
  // finally parse and handle return codes (display help etc...)
  cmd.parseArgs(argc, argv);
  // now query the parsing results

  if (cmd.optionSet("v")) {
    Globals::getInstance()->setVerbose(cmd.getInt("v", 0, 10));
  }
  if (VERB > 0) {
    cerr << extendedGreeter();
  }
  if (cmd.optionSet("Y")) { //NOTE not used?
    //tmpFNs.push_back(cmd.options["Y"]);
    reader->setFile(cmd.options["Y"]);
  }
  if (cmd.optionSet("o")) {
    xmlOutputFN = cmd.options["o"];
  }
  if (cmd.optionSet("Q")) { //NOTE not used?
    parseOptions.calcQuadraticFeatures=true;
  }
  if (cmd.optionSet("y")) {
    Enzyme::setEnzyme(Enzyme::NO_ENZYME);
  }
  if (cmd.optionSet("e")) {
    Enzyme::setEnzyme(Enzyme::ELASTASE);
  }
  if (cmd.optionSet("c")) {
    Enzyme::setEnzyme(Enzyme::CHYMOTRYPSIN);
  }
  if (cmd.optionSet("t")) {
    Enzyme::setEnzyme(Enzyme::THERMOLYSIN);
  }
  if (cmd.optionSet("k")) {
    Enzyme::setEnzyme(Enzyme::PROTEINASEK);
  }
  if (cmd.optionSet("n")) {
    Enzyme::setEnzyme(Enzyme::PEPSIN);
  }
  if (cmd.optionSet("l")) {
    Enzyme::setEnzyme(Enzyme::LYSN);
  }
  if (cmd.optionSet("N")) {
    parseOptions.pngasef=true;
  }
  if (cmd.optionSet("a")) {
    parseOptions.calcAAFrequencies=true;
  }
  if (cmd.optionSet("b")) {
    parseOptions.calcPTMs=true;
  }
  if (cmd.optionSet("m")) {
    int m = cmd.getInt("m", 1, 30000);
    parseOptions.hitsPerSpectrum=m;
  }

  if (cmd.optionSet("2")) {
    spectrumFile = cmd.options["2"];
  }
  if (cmd.optionSet("M")) {
    MassHandler::setMonoisotopicMass(true);
  }
  if (cmd.optionSet("p")) {
    std::vector<std::string> strs;
    boost::split(strs, cmd.options["p"], boost::is_any_of(":,"));
    if (strs.size()<2) {cerr << "Scheme is malformated" << endl; exit(-1);}
    for(unsigned int ix=0; ix+1<strs.size(); ix+=2) {
      parseOptions.ptmScheme[strs[ix][0]]=boost::lexical_cast<int>(strs[ix+1]);   
      if (VERB > 0) {
        cerr << "Interpreting " << strs[ix][0] << " as modification UNIMOD:" << parseOptions.ptmScheme[strs[ix][0]] << endl; 
      }
    }
  }
  
  std::string pattern = "";
  if (cmd.optionSet("P")) pattern = cmd.options["P"];
  bool iscombined = cmd.optionSet("Z");
  
  if(iscombined)
  {
    if(cmd.arguments.size() > 1)
    {
      std::cerr << "Error, there should be only one input file.\n"; 
      exit(-1);
    }
    else if (pattern != "")
    {
      parseOptions.reversedFeaturePattern = pattern;
      parseOptions.iscombined = true;
      targetFN = cmd.arguments[0];
    }
    else
    {
      std::cerr << "Error, pattern should contain a valid set of alphanumberic characters.\n"; 
      exit(-1);
    }
  }
  else
  {
    if (cmd.arguments.size() > 0) targetFN = cmd.arguments[0];
    if (cmd.arguments.size() > 1) decoyFN = cmd.arguments[1];
    
    if(targetFN == "" || decoyFN == "")
    {
      std::cerr << "Error, one of the input files is missing.\n"; 
      exit(-1); 
    }
  }
  // if there are no arguments left...
  if (cmd.arguments.size() == 0) {
      cerr << "Error: too few arguments.";
      cerr << "\nInvoke with -h option for help\n";
      exit(-1); // ...error
  }
  
  
  return true;
}

void msgfdb2Pin::readRetentionTime(string filename) {
  MSReader r;
  Spectrum s;
  r.setFilter(MS2);
  char* cstr = new char[filename.size() + 1];
  strcpy(cstr, filename.c_str());
  // read first spectrum
  r.readFile(cstr, s);
  while(s.getScanNumber() != 0){
    // check whether an EZ lines is available
    if(s.sizeEZ() != 0){
      // for each EZ line (each psm)
      for(int i = 0; i<s.sizeEZ(); i++){
        // save experimental mass and retention time
        scan2rt[s.getScanNumber()].push_back(s.atEZ(i).mh);
        scan2rt[s.getScanNumber()].push_back(s.atEZ(i).pRTime);
      }
    }
    // if no EZ line is available, check for an RTime lines
    else if((double)s.getRTime() != 0){
      scan2rt[s.getScanNumber()].push_back(s.getRTime());
    }
    // if neither EZ nor I lines are available
    else{
      cout << "The ms2 in input does not appear to contain retention time "
          << "information. Please run without -2 option.";
      exit(-1);
    }
    // read next scan
    r.readFile(NULL, s);
  }
  delete[] cstr;
}

void msgfdb2Pin::storeRetentionTime(FragSpectrumScanDatabase* database){
  // for each spectra from the ms2 file
  typedef std::map<int, vector<double> > map_t;
  BOOST_FOREACH(map_t::value_type& i, scan2rt){
    // scan number
    int scanNr = i.first;
    // related retention times
    vector<double>* rTimes = &(i.second);
    if(database->getFSS(scanNr).get()!=0){
      fragSpectrumScan fss = *(database->getFSS(scanNr));
      fragSpectrumScan::peptideSpectrumMatch_sequence& psmSeq =
          fss.peptideSpectrumMatch();
      // retention time to be stored
      double storeMe = 0;
      // if rTimes only contains one element
      if(rTimes->size()==1){
        // take that as retention time
        storeMe = rTimes->at(0);
      }
      else{
        // else, take retention time of psm that has observed mass closest to
        // theoretical mass (smallest massDiff)
        double massDiff = std::numeric_limits<double>::max(); // + infinity
        for (fragSpectrumScan::peptideSpectrumMatch_iterator psmIter_i =
            psmSeq.begin(); psmIter_i != psmSeq.end(); ++psmIter_i) {
          // skip decoy
          if(psmIter_i->isDecoy() != true){
            double cm = psmIter_i->calculatedMassToCharge();
            double em = psmIter_i->experimentalMassToCharge();
            // if a psm with observed mass closer to theoretical mass is found
            if(abs(cm-em) < massDiff){
              // update massDiff
              massDiff = abs(cm-em);
              // get corresponding retention time
              vector<double>::const_iterator r = rTimes->begin();
              for(; r<rTimes->end(); r=r+2){
                double rrr = *r;
                double exm = psmIter_i->experimentalMassToCharge();
                if(*r==psmIter_i->experimentalMassToCharge()){
                  storeMe = *(r+1);
                  r = rTimes->end();
                }
              }
            }
          }
        }
      }
      // store retention time for all psms in fss
      for (fragSpectrumScan::peptideSpectrumMatch_iterator psmIter =
          psmSeq.begin(); psmIter != psmSeq.end(); ++psmIter) {
        psmIter->observedTime().set(storeMe);
      }
      database->putFSS(fss);
    }
  }
}

int msgfdb2Pin::run() {
  // Content of msgfdb files is merged: preparing to write it to xml file
  ofstream xmlOutputStream;
  xmlOutputStream.open(xmlOutputFN.c_str());
  if(!xmlOutputStream && xmlOutputFN != ""){
    cerr << "ERROR: invalid path to output file: " << xmlOutputFN << endl;
    cerr << "Please invoke msgfdb2pin with a valid -o option" << endl;
    exit(-1);
  }
  xercesc::XMLPlatformUtils::Initialize ();

  // initializing xercesc objects corresponding to pin element...
  // ... <featureDescriptions>
  std::auto_ptr<percolatorInNs::featureDescriptions>
  fdes_p (new ::percolatorInNs::featureDescriptions());

  // ... <process_info>
  percolatorInNs::process_info::command_line_type command_line = call;
  std::auto_ptr<percolatorInNs::process_info>
  proc_info (new ::percolatorInNs::process_info(command_line));

  // ... <experiment>
  std::auto_ptr< ::percolatorInNs::experiment >
  ex_p (new ::percolatorInNs::experiment("mitt enzym", proc_info, fdes_p));

  int maxCharge = -1;
  int minCharge = 10000;

  vector<FragSpectrumScanDatabase*> databases;

  if (!parseOptions.iscombined) {
    // First we only search for the maxCharge and minCharge.
    // This done by passing the argument justSearchMaxMinCharge
    
    std::cerr << "Reading input from msgfdb files:\n";
    
    // Search and set max and min charge
    reader->translateFileToXML(targetFN,ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), false /* is_decoy */, parseOptions,
        &maxCharge, &minCharge,  Reader::justSearchMaxMinCharge, databases,
        0/*, tmpDirs, tmpFNs*/);
    reader->translateFileToXML(decoyFN, ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), true /* is_decoy */, parseOptions,
        &maxCharge, &minCharge, Reader::justSearchMaxMinCharge, databases,
        0/*, tmpDirs, tmpFNs*/);
    // Now we do full parsing of the msgfdb file, and translating it to XML
    reader->translateFileToXML(targetFN,ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), false /* is_decoy */, parseOptions,
        &maxCharge, &minCharge,  Reader::fullParsing, databases,
        0/*, tmpDirs, tmpFNs*/);
    reader->translateFileToXML(decoyFN, ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), true /* is_decoy */, parseOptions,
        &maxCharge, &minCharge, Reader::fullParsing, databases,
        0/*, tmpDirs, tmpFNs*/);

  } else {
    // First we only search for the maxCharge and minCharge.
    //This done by passing the argument justSearchMaxMinCharge
    
    std::cerr << "Reading input from a combined (target-decoy) msgfdb file .." << std::endl;
    
     // Search and set max and min charge
    reader->translateFileToXML(targetFN,ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), false /* is_decoy */, parseOptions,
        &maxCharge, &minCharge, Reader::justSearchMaxMinCharge, databases,
        0/*, tmpDirs, tmpFNs*/);
    // Now we do full parsing of the msgfdb file, and translating it to XML
    reader->translateFileToXML(targetFN,ex_p->featureDescriptions(),
        ex_p->fragSpectrumScan(), false /* is_decoy */, parseOptions,
        &maxCharge, &minCharge, Reader::fullParsing, databases,
        0/*, tmpDirs, tmpFNs*/);
  }

  // read retention time if msgfdb2pin was invoked with -2 option
  if (spectrumFile.size() > 0) {
    readRetentionTime(spectrumFile);
    databases[0]->initRTime(&scan2rt);
    storeRetentionTime(databases[0]);
  }

  xercesc::XMLPlatformUtils::Terminate();

  string schema_major = boost::lexical_cast<string>(PIN_VERSION_MAJOR);
  string schema_minor = boost::lexical_cast<string>(PIN_VERSION_MINOR);
  string headerStr = "<?xml version=\"1.0\" encoding=\"UTF-8\"?> \n" +
      string("<experiment xmlns=\"") + PERCOLATOR_IN_NAMESPACE + "\"" +
      " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"" +
      " xsi:schemaLocation=\"" + PERCOLATOR_IN_NAMESPACE +
      " https://github.com/percolator/percolator/raw/pin-" + schema_major +
      "-" + schema_minor + "/src/xml/percolator_in.xsd\"> \n";
  if (xmlOutputFN == "") cout << headerStr;
  else {
    xmlOutputStream << headerStr;
    cerr <<  "The output will be written to " << xmlOutputFN << endl;
  }

  string enzymeStr = "\n<enzyme>" + Enzyme::getStringEnzyme() + "</enzyme>\n";
  if (xmlOutputFN == "") cout << enzymeStr;
  else xmlOutputStream << enzymeStr;

  string commandLine = "\n<process_info>\n" +
      string("  <command_line>") + call.substr(0,call.length()-1)
      + "</command_line>\n" + "</process_info>\n";
  if (xmlOutputFN == "") cout << commandLine;
  else xmlOutputStream << commandLine;

  xercesc::XMLPlatformUtils::Initialize ();

  cerr << "\nWriting output:\n";
  // print to cout (or populate xml file)
  // print features
  {
    serializer ser;
    if (xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    ser.next ( PERCOLATOR_IN_NAMESPACE, "featureDescriptions",
        ex_p->featureDescriptions());
  }

  // print fragSpecturmScans
  std::cerr << "Databases : " << databases.size() << std::endl;
  for(int i=0; i<databases.size();i++) {
    serializer ser;
    if (xmlOutputFN == "") ser.start (std::cout);
    else ser.start (xmlOutputStream);
    if(VERB>1){
      cerr << "outputting content of " << databases[i]->id
          << " (and correspondent decoy file)\n";
    }
    databases[i]->print(ser);
    databases[i]->terminte();
  }

  // print closing tag
  if (xmlOutputFN == "") std::cout << "</experiment>" << std::endl;
  else {
    xmlOutputStream << "</experiment>" << std::endl;
    xmlOutputStream.close();
  }

  xercesc::XMLPlatformUtils::Terminate();

  cerr << "\nAll the input files have been successfully processed"<< endl;

  return 0;
}

int main(int argc, char** argv) {
  msgfdb2Pin* pmsgfdb2Pin = new msgfdb2Pin();
  int retVal = -1;
  if (pmsgfdb2Pin->parseOpt(argc, argv)) {
    retVal = pmsgfdb2Pin->run();
  }
  delete pmsgfdb2Pin;
  Globals::clean();
  return retVal;
}