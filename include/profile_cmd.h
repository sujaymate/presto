#ifndef __profile_cmd__
#define __profile_cmd__
/*****
  command line parser interface 
                        -- generated by clig (Version: 1.1.3)

  The command line parser `clig':
  (C) 1995 Harald Kirsch (kir@iitb.fhg.de)
*****/

#ifndef FALSE
#  define FALSE (1==0)
#  define TRUE  (!FALSE)
#endif

typedef struct s_Cmdline {
  /***** -phs: Offset phase for the profile */
  char phsP;
  double phs;
  int phsC;
  /***** -p: The folding period (s) */
  char pP;
  double p;
  int pC;
  /***** -pd: The folding period derivative (s/s) */
  char pdP;
  double pd;
  int pdC;
  /***** -pdd: The folding period 2nd derivative (s/s^2) */
  char pddP;
  double pdd;
  int pddC;
  /***** -f: The folding frequency (hz) */
  char fP;
  double f;
  int fC;
  /***** -fd: The folding frequency derivative (hz/s) */
  char fdP;
  double fd;
  int fdC;
  /***** -fdd: The folding frequency 2nd derivative (hz/s^2) */
  char fddP;
  double fdd;
  int fddC;
  /***** -n: The number of bins in the profile.  Defaults to the number of sampling bins which correspond to one folded period */
  char proflenP;
  int proflen;
  int proflenC;
  /***** -psr: Name of pulsar to fold (do not include J or B) */
  char psrnameP;
  char* psrname;
  int psrnameC;
  /***** -rzwcand: The candidate number to fold from 'infile'_rzw.cand */
  char rzwcandP;
  int rzwcand;
  int rzwcandC;
  /***** -rzwfile: Name of the rzw search file to use (include the full name of the file) */
  char rzwfileP;
  char* rzwfile;
  int rzwfileC;
  /***** -bincand: Fold a binary pulsar but take the input data from this candidate number in 'infile'_bin.cand */
  char bincandP;
  int bincand;
  int bincandC;
  /***** -onoff: A list of white-space separated pairs of numbers from 0.0 to 1.0 that designate barycentric times in our data set when we will actually keep the data. (i.e. '-onoff 0.1 0.4 0.7 0.9' means that we will fold the data set during the barycentric times 0.1-0.4 and 0.7-0.9 of the total time length of the data set) */
  char onoffP;
  char* onoff;
  int onoffC;
  /***** -bin: Fold a binary pulsar.  Must include all of the following parameters */
  char binaryP;
  /***** -pb: The orbital period (s) */
  char pbP;
  double pb;
  int pbC;
  /***** -x: The projected orbital semi-major axis (lt-sec) */
  char asinicP;
  double asinic;
  int asinicC;
  /***** -e: The orbital eccentricity */
  char eP;
  double e;
  int eC;
  /***** -To: The time of periastron passage (MJD) */
  char ToP;
  double To;
  int ToC;
  /***** -w: Longitude of periastron (deg) */
  char wP;
  double w;
  int wC;
  /***** -wdot: Rate of advance of periastron (deg/yr) */
  char wdotP;
  double wdot;
  int wdotC;
  /***** -xwin: Send graphics output to the screen */
  char xwinP;
  /***** -ps: Send graphics output to a Postscript file */
  char psP;
  /***** -both: Send graphics output both the screen and a Postscript file */
  char bothP;
  /***** -disp: Don't calculate a new profile.  Just display a previously calculated profile in 'infile'.prof.  Must be called with either -ps or -xwin */
  char dispP;
  /***** -mak: Determine folding parameters from 'infile.mak' */
  char makefileP;
  /***** -noerr: Do not plot error bars */
  char noerrP;
  /***** uninterpreted command line parameters */
  int argc;
  /*@null*/char **argv;
  /***** the whole command line concatenated */
  char *full_cmd_line;
} Cmdline;


extern char *Program;
extern void usage(void);
extern /*@shared*/Cmdline *parseCmdline(int argc, char **argv);

extern void showOptionValues(void);

#endif

