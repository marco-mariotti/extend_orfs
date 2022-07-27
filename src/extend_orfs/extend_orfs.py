#! /usr/bin/env python2.7
from string import *
import sys
from commands import *
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *

help_msg="""Program to extend protein coding predictions upstream and/or downstream. The aim is to include the closest stop codon downstream; and the most distant methionine upstream before encountering a stop. Some genes in output may still be missing the stop or start: for the 5' extension, this happens if a stop codon is encountered before any start; for the stop, this may happen if the gene prediction reaches the very end of a sequence entry. 

extend_orfs.py  -i genes.gff  -t target_file.fa  -o output_file.gff  [options]    

### Options:
-pep     +      peptide sequence output file in fasta format  
-cds     +      coding sequence output file in fasta format  
-tag     +      normally only the "CDS" lines are read from the input file. This is to specify another tag. You can use "*" to accept any tag.
-starts  +      provides comma separated list of codons accepted as starts, in DNA format. Defaults to ATG only. Example: -starts ATG,TTG
-stops   +      provides comma separated list of codons accepted as stops, in DNA format. Defaults to TAG,TAA,TGA
-only_down      extends only downstream, not upstream   
-only_up        extends only upstream, not downstream   
-genetic_code + NCBI id of genetic code used for translation. This also changes the stop codons, unless they are provided on the command line
-print_opt      print currently active options
-h OR --help    print this help and exit"""

command_line_synonyms={}

def_opt= { #'temp':'/users-d3/mmariotti/temp', 
'i':0, 't':0, 'o':0,
'tag':'cds', 'starts':0, 'stops':0,
'genetic_code':0,
'pep':0, 'cds':0,
'only_up':0, 'only_down':0,
'v':0,
}


#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 

  if not opt['o']: raise Exception, "ERROR you must choose an output file with option -o"
  ## checking input
  input_file=opt['i'];       check_file_presence(input_file, 'input_file')
  tag=opt['tag']
  printerr('Loading genes with tag "{0}" in file: {1} ...'.format(tag, input_file), 1)
  all_genes= load_all_genes(input_file, tag=tag, keep_program=True)
  printerr('-- Loaded {0} gene structures.'.format(len(all_genes)), 1)

  target_file=opt['t'];       check_file_presence(target_file, 'target_file')
  printerr('Loading sequence database in memory: {0} ...'.format(target_file), 1)
  load_sequence_db(target_file)
  printerr('-- Loaded {0} sequences.'.format(len(get_MMlib_var('sequence_db').keys())), 1)

  keep_sequence_in_gene_object=  opt['cds'] or opt['pep']

  if opt['genetic_code']:
    set_genetic_code(opt['genetic_code'])

  if not opt['starts'] : opt['starts'] ='ATG'
  if not opt['stops']  : 
    if opt['genetic_code']:      
      opt['stops']=','.join( [c for c,aa in  get_genetic_code_table().iteritems() if aa=='*'] )
    else: opt['stops']  ='TAG,TAA,TGA'  

  allowed_stops={}; allowed_starts={}
  for stop in  opt['stops'].split(','):  allowed_stops[stop]=1
  for start in opt['starts'].split(','): allowed_starts[start]=1


  extend_up= not opt['only_down'];   extend_down= not opt['only_up']
  how_many_extended=0
  outfile_h=open(opt['o'], 'w')
  if opt['cds']:  cds_outfile_h=open(opt['cds'], 'w')
  if opt['pep']:  pep_outfile_h=open(opt['pep'], 'w')
  for g in all_genes:
    if not g.chromosome in get_MMlib_var('sequence_db'): raise Exception, "ERROR chromosome/scaffold {0} not found in target file {1}".format(g.chromosome, target_file)
    service('extending orf of gene {0} ... '.format(g.id))
    bounds_before=g.boundaries()
    g.extend_orf(chromosome_length=len(get_MMlib_var('sequence_db')[g.chromosome]), stops=allowed_stops, starts=allowed_starts, up=extend_up, down=extend_down, get_seq=lambda x:replace(upper(x.fast_sequence()),'U','T'),  extension_parameter=1000, keep_seq=keep_sequence_in_gene_object)
    bounds_after= g.boundaries()
    if bounds_before != bounds_after: how_many_extended+=1
    print >> outfile_h, g.gff(tag=tag)
    if opt['cds']: print >> cds_outfile_h, ">{0}\n{1}".format(g.id, fasta(g.seq))
    if opt['pep']: print >> pep_outfile_h, ">{0}\n{1}".format(g.id, fasta(transl(g.seq)))
  
  printerr('-- Finished. {0} gene structures processed, {1} were extended.'.format(len(all_genes), how_many_extended) , 1)  
  ###############



#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass

  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
