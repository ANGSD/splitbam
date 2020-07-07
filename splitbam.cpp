//gpl thorfinn@binf.ku.dk
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <getopt.h>
#include <cassert>
#include <ctime>

char out_mode[5]="wb";

htsFormat *dingding2 =(htsFormat*) calloc(1,sizeof(htsFormat));

int usage(FILE *fp, int is_long_help)
{
    fprintf(fp,
"\n"
"Usage: ./splitbam [options] <in.bam>|<in.sam>|<in.cram> \n"
"\n"
"Options:\n"
// output options
"  -b       output BAM\n"
"  -C       output CRAM (requires -T)\n"
"  -o FILE  output file name \n"
"  -T FILE  reference in the fastaformat (required from reading and writing crams)\n"
"  -@ INT   Number of threads to use\n"
"  -q INT   Mapping quality filter\n"
"  -m       Discard unmapped reads (default off)\n"
"  -v       Verbose mode\n"
	    );
    fprintf(fp,
	    "\nNotes:\n"
	    "\n"
	    "1. This program is usefull for splitting a sorted bam/cram into multiple files\n");

    return 0;
}

void parse_sequencingdata(char *fn_out,char *refName,char *fname,int nthreads,int mapped_only,int mapq,int ineach){
  htsThreadPool p = {NULL, 0};
  samFile *in=NULL;
  samFile *out=NULL;
  char onam1[2048]="";
  int whichfile=0;
  char lastnamewritten[2048];
  if(refName){
    char *ref =(char*) malloc(10 + strlen(refName) + 1);
    sprintf(ref, "reference=%s", refName);
    hts_opt_add((hts_opt **)&dingding2->specific,ref);
    free(ref);
  }
  
  if((in=sam_open_format(fname,"r",dingding2))==NULL ){
    fprintf(stderr,"[%s] nonexistant file: %s\n",__FUNCTION__,fname);
    exit(0);
  }

  bam_hdr_t  *hdr = sam_hdr_read(in);

  if(strstr(fname,".cram")!=NULL &&out_mode[1]=='c'&&refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    exit(0);
  }
  
  if(out_mode[1]=='c'&&refName==NULL){
    fprintf(stderr,"\t-> cram file requires reference with -T FILE.fa \n");
    exit(0);
  }

 
  if(nthreads>1){
    if (!(p.pool = hts_tpool_init(nthreads))) {
      fprintf(stderr, "Error creating thread pool\n");
      exit(0);
    }
    hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
  }
    

  
  bam1_t *b = bam_init1();

  int ret;

  int nproc;
  int nwritten =0;
  while(((ret=sam_read1(in,hdr,b)))>0){
    if(nwritten>ineach){
      if(strcmp(lastnamewritten,bam_get_qname(b))!=0){
	 assert(sam_close(out)==0);
	 out=NULL;
	 nwritten=0;
      }
    }
    nproc++;
    if(out==NULL){
      nwritten = 0;
      if(out_mode[1]=='b')
	snprintf(onam1,2024,"%s.%d.bam",fn_out,whichfile++);
      else
	snprintf(onam1,2024,"%s.%d.cram",fn_out,whichfile++);

      if ((out = sam_open_format(onam1, out_mode, dingding2)) == 0) {
	fprintf(stderr,"Error opening file for writing\n");
	exit(0);
      }
      if (nthreads>1&&out) hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);
      assert(sam_hdr_write(out, hdr) == 0);
    }
   
    if(mapped_only!=0){
      if(b->core.flag&4)
	continue;
    }
    		    
    if(mapq!=-1 && b->core.qual<mapq)
      continue;
    strcpy(lastnamewritten,bam_get_qname(b));
    nwritten++;
    assert(sam_write1(out, hdr,b)); //write into the file containing the pcrdups+normal reads
    
  }
  assert(sam_close(out)==0);
  assert(sam_close(in)==0);
  
  bam_destroy1(b);
  hts_opt_free((hts_opt *)dingding2->specific);
  free(dingding2);
}


int main(int argc, char **argv){
  //  int VERBOSE = 0;

  clock_t t=clock();
  time_t t2=time(NULL);

  char *fname,*refName;

  
  fname=refName=NULL;
  char *fn_out = NULL;
  int c;
  int nthreads = 1;

  int mapq =-1;
  int mapped_only = 0;
  int ineach = 5000000;//5mio reads in one file
  if(argc==1){
    usage(stdout,0);
    return 0;
  }
  //fix these
  static struct option lopts[] = {
    {"add", 1, 0, 0},
    {"append", 0, 0, 0},
    {"delete", 1, 0, 0},
    {"verbose", 0, 0, 0},
    {"create", 1, 0, 'c'},
    {"file", 1, 0, 0},
    {NULL, 0, NULL, 0}
  };
  
  while ((c = getopt_long(argc, argv,
			  "bCo:T:@:q:m:n:v",
			  lopts, NULL)) >= 0) {
    switch (c) {
    case 'b': out_mode[1] = 'b'; break;
    case 'C': out_mode[1] = 'c'; break;
    case 'T': refName = strdup(optarg); break;
    case 'o': fn_out = strdup(optarg); break;
    case '@': nthreads = atoi(optarg); break;
    case 'q': mapq = atoi(optarg); break;
    case 'n': ineach = atoi(optarg); break;
    case 'm': mapped_only = 1; break;
      //    case 'v': VERBOSE = 1; break;
        case '?':
	  if (optopt == '?') {  // '-?' appeared on command line
	    return usage(stdout,0);
	  } else {
	    if (optopt) { // Bad short option
	      fprintf(stdout,"./superduper invalid option -- '%c'\n", optopt);
	    } else { // Bad long option
	      // Do our best.  There is no good solution to finding
	      // out what the bad option was.
	      // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
	      if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
		fprintf(stdout,"./superduper unrecognised option '%s'\n",argv[optind - 1]);
	      }
	    }
	    return 0;//usage(stderr, 0);
	  }
    default:
      fprintf(stderr,"adsadsfasdf\n");
      fname = strdup(optarg);
      fprintf(stderr,"assinging: %s to fname:%s\n",optarg,fname);
      break;
    }
  }
  if(optind<argc)
    fname = strdup(argv[optind]);
  
  if(!fname){
    fprintf(stderr,"\t-> No input file specified\n");
    usage(stdout,0);
    return 0;
  }

  if(!fn_out){
    fprintf(stderr,"\t-> No output file specified\n");
    usage(stdout,0);
    return 0;
  }
  
  if(fname)
    parse_sequencingdata(fn_out,refName,fname,nthreads,mapped_only,mapq,ineach);
 
  fprintf(stderr,
	  "\t[ALL done] cpu-time used =  %.2f sec\n"
	  "\t[ALL done] walltime used =  %.2f sec\n"
	  ,(float)(clock() - t) / CLOCKS_PER_SEC, (float)(time(NULL) - t2));  
  return 0;
}

