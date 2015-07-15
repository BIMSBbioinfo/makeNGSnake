import os

if not os.path.exists("./output/FASTQ"):
    os.makedirs("./output/FASTQ")
if not os.path.exists("./output/ENA_reports"):
    os.makedirs("./output/ENA_reports")
     
     
# In case this is being run on a cluster then make sure that it is executed on the host node
localrules: download

EXP_IDs = ["ERX620541"] #can be also study accession, run accession etc.
LIBRARY_STRATEGY = "ChIP-Seq"


def wget_file(file_path, file_url, filename="the requested file"):
    """
    Download a file named 'filename' from a url 'file_url' to 'file_path'  
    """
    if(os.path.isfile(file_path)):
        print('Using locally cached version of '+ filename +' found here:\n' + file_path)
    else:
        shell("wget -O " + file_path + " " + file_url)
        
def url_ENA(id):
    return "http://www.ebi.ac.uk/ena/data/warehouse/filereport\?accession\="+id+"\&result=read_run\&fields=study_accession,secondary_study_accession,sample_accession,secondary_sample_accession,experiment_accession,run_accession,tax_id,scientific_name,instrument_model,library_layout,library_strategy,library_selection,fastq_ftp"

                       
def links_runs(exp_ids, library_strategy):
    """function that extracts ALL links to fastq files of runs from ALL experiments and
       returns a dict with key=name_of_file and value = link
    """    
    dct = {}
    for exp_id in exp_ids:
        #download file that contrains info about url to fastq files of runs within experiment exp_id
        wget_file("output/ENA_reports/"+exp_id+".txt", url_ENA(exp_id), exp_id+".txt")
        #extract urls of fastq file of runs and put them into dictionary dct
        f = open("output/ENA_reports/"+exp_id+".txt", 'r')
        l = 0 
        header = []
        for line in f:
              line = line.rstrip().split("\t")
              if l==0:
                  header=line
                  l+=1
                  continue 
              d = dict(zip(header, line))
              if d["library_strategy"]==library_strategy:
                      fastaq_paths = d["fastq_ftp"].split(";")
                      for i in fastaq_paths:
                         runid = i.split("/")[-1].split(".")[0]
                         dct[runid] = i
    if not bool(dct): #if dict is empty
      print("There are no runs..")
    return dct




RUN_ID_URL = links_runs(EXP_IDs, LIBRARY_STRATEGY)



rule final: 
     input: 
         expand("output/FASTQ/{run_id}.fastq.gz", run_id=RUN_ID_URL.keys())


# Get fastq files from the ENA database.
rule download: 
     output: 
         "output/FASTQ/{run_id}.fastq.gz" 
     run: 
          wget_file("./output/FASTQ/"+wildcards.run_id+".fastq.gz", RUN_ID_URL[wildcards.run_id], wildcards.run_id+".fastq.gz")
 
  
rule clean:
    """Clean output directory. It will remove all files from output/FASTQ and output/ENA_reports directories
    """
    shell: "rm -r output/FASTQ/*; rm -r output/ENA_reports/*"                              
 
 
 
                            