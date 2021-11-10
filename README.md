# rattleFilter

When using RATTLE (https://github.com/comprna/RATTLE) it might be possible to obtain multi overlapping small pieces of transcripts.
Use this script to keep the best transcripts and filter the transcriptome obtained with RATTLE.




## Usage

rattleFilter.py [-h] [-g RATTLE_GTF] [-f RATTLE_FASTQ] [-o OUTPUT]
                     [-n MINNUCLEOTIDES] [-i INTERSECTLENGTH]
                     [-p PERCENTAGEREADS] [-m MINREADLENGTH]
                     [-M MAXREADLENGTH]

optional arguments:
  -h, --help           
  
  show this help message and exit

  -g RATTLE_GTF, --RATTLE_GTF RATTLE_GTF
  
                        RATTLE transcripts mapped to reference genome in GTF
                        format

  -f RATTLE_FASTQ, --RATTLE_FASTQ RATTLE_FASTQ
  
                        RATTLE output 'transcriptome.fq' in FASTQ format

  -o OUTPUT, --output OUTPUT
  
                        specify an output directory (default: current folder)

  -n MINNUCLEOTIDES, --minNucleotides MINNUCLEOTIDES
  
                        cutoff of transcript length (nt) (default: 100)

  -i INTERSECTLENGTH, --intersectLength INTERSECTLENGTH
  
                        minimum percentage of transcript overlap to be
                        accepted in a group of transcripts (default: 0.5)

  -p PERCENTAGEREADS, --percentageReads PERCENTAGEREADS
  
                        minimum percentage of associated reads that a
                        transcript (in a group of overlapping transcripts)
                        must have if it is not the longest one (default: 0.5)

  -m MINREADLENGTH, --minReadLength MINREADLENGTH
  
                        cutoff of MINUMUM length that transcript selected by
                        minumum percentage of associated reads (-p) should
                        have. It is a length percentage of the transcript with
                        the highest number of associated reads (default: 0.5)

  -M MAXREADLENGTH, --maxReadLength MAXREADLENGTH
  
                        cutoff of MAXIMUM length that transcript selected by
                        minumum percentage (-p) should have. It is a length
                        percentage of the transcript with the highest number
                        of associated reads (default: 2)
                        
                        
## Requirements
python version 3             
                    
