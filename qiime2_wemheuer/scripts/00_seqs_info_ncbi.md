# Поиск данных секвенирования ассоциированных с проектом PRJNA594470 через CLI
# Aмпликонное секвенирование с универасальными праймерами 16S, ITS1 (отличные от моих)
# Эндофитные микроорганизмы Theobroma cacao через


# Поиск данных секвенирования ассоциированных с проектом PRJNA925518 через CLI
# Ампликонное секвенирование эндофитов из листьев какао устойчивого и чувствительного к фитофторе сортов


# Основная таблица всех ранов — размеры, имена, платформа, стратегия
```bash
esearch -db sra -query "PRJNA594470" | \
    efetch -format runinfo > PRJNA594470_runinfo.csv
```
git
# Посмотреть: Run, size_MB, LibraryStrategy, LibraryLayout, spots, SampleName и т.д.
```bash
column -t -s',' PRJNA594470_runinfo.csv | less -S
```
qi  
# Сколько ранов всего
```bash
esearch -db sra -query "PRJNA594470" | efetch -format runinfo | tail -n +2 | wc -l
# 128
```

# Суммарный размер в GB
```bash
awk -F',' 'NR>1 {sum+=$8} END {print sum/1024 " GB"}' PRJNA594470_runinfo.csv
# 3.71777 GB
```

Я утсановил pysrddb чтобы посмотреть метаданные. 
```bash
pip install pysradb   # если нет
```

# Все метаданные с аттрибутами образцов
```bash
pysradb metadata --detailed PRJNA594470 > PRJNA594470_metadata.tsv
```
# Только акцессии ранов
```bash
pysradb srp-to-srr PRJNA594470
```

# Извлечь только SRR-акцессии
```bash
cut -d',' -f1 PRJNA594470_runinfo.csv | grep "^SRR" > runs.txt
wc -l runs.txt
# 128 runs.txt
```

## РАЗВЕДКА ##
```bash
pwd
# /mnt/c/Users/nmmak/Desktop/cacao_flower_microbiome
mamba activate bioinf

mkdir data_wemheuer
esearch -db sra -query "PRJNA594470" | efetch -format runinfo >  data_wemheuer/runinfo.csv

head -3 data_wemheuer/runinfo.csv
# Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
# SRR10614063,2020-03-12 09:40:51,2019-12-10 04:53:41,73219,43605277,73219,595,29,,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614063.sralite.1,SRX7293343,Fun_T32,AMPLICON,PCR,METAGENOMIC,PAIRED,0,0,ILLUMINA,Illumina MiSeq,SRP235361,PRJNA594470,,594470,SRS5786637,SAMN13519336,simple,1297858,leaf metagenome,Tree32,,,,,,,no,,,,,GENOMIC AND APPLIED MICROBIOLOGY & GOETTINGEN GENOMICS LABORATORY,SRA1008811,,public,3678FF6290A9F27D1C3FBB95F3B0D8EF,0EC9187FFA58135EF4E869DBE14AAB35
# SRR10614062,2020-03-12 09:40:51,2019-12-10 04:53:16,123636,74168797,123636,599,46,,https://sra-downloadb.be-md. ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614062.sralite.1,SRX7293344,Fun_T33,AMPLICON,PCR,METAGENOMIC,PAIRED,0,0,ILLUMINA,Illumina MiSeq,SRP235361,PRJNA594470,,594470,SRS5786638,SAMN13519337,simple,1297858,leaf metagenome,Tree33,,,,,,,no,,,,,GENOMIC AND APPLIED MICROBIOLOGY & GOETTINGEN GENOMICS LABORATORY,SRA1008811,,public,CDA6478EC7406E77FA31F5169A225B40,C32BD2E8DEA4869D9F76BD5AE690AE14

wc -l data_wemheuer/runinfo.csv 
# 129 runinfo.csv

awk -F',' 'NR>1 {sum+=$8} END {printf "Total: %.1f GB\n", sum/1024}' data_wemheuer/runinfo.csv
# Total: 3.7 GB

pysradb metadata --detailed PRJNA594470 > data_wemheuer/metadata.tsv
head -3 data_wemheuer/metadata.tsv 
# run_accession   study_accession study_title     experiment_accession    experiment_title        experiment_desc organism_taxid  organism_name   library_name    library_strategy        library_source  library_selection       library_layoutsample_accession        sample_title    biosample       bioproject      instrument      instrument_model        instrument_model_desc   total_spots     total_size      run_total_spots run_total_bases run_alias       public_filename public_size   public_date     public_md5      public_version  public_semantic_name    public_supertype        public_sratoolkit       aws_url aws_free_egress aws_access_type public_url      ncbi_url        ncbi_free_egress        ncbi_access_type      gcp_url gcp_free_egress gcp_access_type experiment_alias        collection_date env_broad_scale env_local_scale env_medium      geo_loc_name    host    lat_lon altitude        temp    tree_height     plantation      replicate    biosamplemodel   ena_fastq_http  ena_fastq_http_1        ena_fastq_http_2        ena_fastq_ftp   ena_fastq_ftp_1 ena_fastq_ftp_2
# SRR10614059     SRP235361       Theobroma cacao leaf metagenome SRX7293347      Fungal endophyte communities in leaves of Theobroma cacao, sample Tree36        Fungal endophyte communities in leaves of Theobroma cacao, sample Tree36     1297858  leaf metagenome Fun_T36 AMPLICON        METAGENOMIC     PCR     PAIRED  SRS5786641      -       SAMN13519340    PRJNA594470     Illumina MiSeq  Illumina MiSeq  ILLUMINA        76746   29965293        76746   45406069        Fun_T36_R1.fastq.gz   SRR10614059.sralite     11694441        2020-06-24 01:50:39     798c2a53ffd3a5a32dcabcf015ce2786        1       SRA Lite        Primary ETL     1       s3://sra-pub-zq-8/SRR10614059/SRR10614059.sralite.1     s3.us-east-1 aws identity     https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR010/614/SRR10614059.sralite.1        https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR010/614/SRR10614059.sralite.1        worldwide    anonymous        gs://sra-pub-zq-108/SRR10614059/SRR10614059.zq.1        gs.us-east1     gcp identity    -       2014    cropland biome  cacao plantation        leaf    Cameroon: Bakoa Theobroma cacao 4.5739 N 11.1798 E      683m    24.878.3m     Bakoa   Bakoa_rep4      MIGS/MIMS/MIMARKS.plant-associated      -       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/059/SRR10614059/SRR10614059_1.fastq.gz       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/059/SRR10614059/SRR10614059_2.fastq.gz     -       era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/059/SRR10614059/SRR10614059_1.fastq.gz    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/059/SRR10614059/SRR10614059_2.fastq.gz
# SRR10614060     SRP235361       Theobroma cacao leaf metagenome SRX7293346      Fungal endophyte communities in leaves of Theobroma cacao, sample Tree35        Fungal endophyte communities in leaves of Theobroma cacao, sample Tree35     1297858  leaf metagenome Fun_T35 AMPLICON        METAGENOMIC     PCR     PAIRED  SRS5786640      -       SAMN13519339    PRJNA594470     Illumina MiSeq  Illumina MiSeq  ILLUMINA        10164   4085271 10164   5967227 Fun_T35_R1.fastq.gz  SRR10614060.sralite      1606168 2020-06-20 10:57:53     a9fb63a3500fe50a5e6fa259d7b54c86        1       SRA Lite        Primary ETL     1       s3://sra-pub-zq-7/SRR10614060/SRR10614060.sralite.1     s3.us-east-1    aws identity    https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614060.sralite.1     https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614060.sralite.1       worldwide       anonymous       gs://sra-pub-zq-108/SRR10614060/SRR10614060.zq.1      gs.us-east1     gcp identity    -       2014    cropland biome  cacao plantation        leaf    Cameroon: Bakoa Theobroma cacao 4.5739 N 11.1798 E      683m    24.87   7.8m    Bakoa   Bakoa_rep3    MIGS/MIMS/MIMARKS.plant-associated      -       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/060/SRR10614060/SRR10614060_1.fastq.gz       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/060/SRR10614060/SRR10614060_2.fastq.gz       -    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/060/SRR10614060/SRR10614060_1.fastq.gz     era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/060/SRR10614060/SRR10614060_2.fastq.gz
```

```bash
mkdir data_schmidt
esearch -db sra -query "PRJNA925518" | efetch -format runinfo >  data_schmidt/runinfo.csv

head -3 data_schmidt/runinfo.csv
# Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
# SRR23126323,2023-09-07 15:36:33,2023-01-19 15:03:17,169266,75999907,169266,448,16,,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR023/23126/SRR23126323/SRR23126323.lite.1,SRX19077291,m1_1_9_16S,AMPLICON,PCR,METAGENOMIC,PAIRED,0,0,ILLUMINA,Illumina NovaSeq 6000,SRP418302,PRJNA925518,,925518,SRS16492251,SAMN32798684,simple,662107,phyllosphere metagenome,m1_1_9,,,,,,,no,,,,,MARS WRIGLEY,SRA1576723,,public,676D103FD6001CA8017824E5A2C8FD9E,6815A2E04760F05A17183B15941CF31B
# SRR23126322,2023-09-07 15:36:33,2023-01-19 15:03:24,171398,76957320,171398,448,15,,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR023/23126/SRR23126322/SRR23126322.lite.1,SRX19077292,m1_2_1_16S,AMPLICON,PCR,METAGENOMIC,PAIRED,0,0,ILLUMINA,Illumina NovaSeq 6000,SRP418302,PRJNA925518,,925518,SRS16492252,SAMN32798685,simple,662107,phyllosphere metagenome,m1_2_1,,,,,,,no,,,,,MARS WRIGLEY,SRA1576723,,public,0F503B3799848E593E03E9DFB3EC3EC1,F4F069A5596A8DD2183AD9D75646EFA7

wc -l data_schmidt/runinfo.csv 
# 80 runinfo.csv

awk -F',' 'NR>1 {sum+=$8} END {printf "Total: %.1f GB\n", sum/1024}' data_schmidt/runinfo.csv
# Total: 1.5 GB

pysradb metadata --detailed PRJNA925518 > data_schmidt/metadata.tsv
head -3 data_schmidt/metadata.tsv
# run_accession   study_accession study_title     experiment_accession    experiment_title        experiment_desc organism_taxid  organism_name   library_name    library_strategy        library_source  library_selection       library_layout  sample_accession        sample_title    biosample       bioproject      instrument      instrument_model        instrument_model_desc   total_spots     total_size      run_total_spots run_total_bases run_alias       public_filename public_size     public_date     public_md5      public_version  public_semantic_name    public_supertype        public_sratoolkit       aws_url aws_free_egress aws_access_typepublic_url       ncbi_url        ncbi_free_egress        ncbi_access_type        gcp_url gcp_free_egress gcp_access_type experiment_alias        collection_date env_broad_scale env_local_scale env_medium      geo_loc_name    host    lat_lon id      biosamplemodel  ena_fastq_http  ena_fastq_http_1        ena_fastq_http_2        ena_fastq_ftp   ena_fastq_ftp_1 ena_fastq_ftp_2
# SRR23126322     SRP418302       Cacao phyllosphere microbiome is related to Phytophthora palmivora resistance   SRX19077292     amplicon sequencing of Theobroma cacao: phyllosphere microbiome amplicon sequencing of Theobroma cacao: phyllosphere microbiome 662107  phyllosphere metagenome m1_2_1_16S      AMPLICON        METAGENOMIC     PCR     PAIRED  SRS16492252     -       SAMN32798685    PRJNA925518     Illumina NovaSeq 6000  Illumina NovaSeq 6000    ILLUMINA        171398  16358193        171398  76957320        m1_2_1_1.fq     SRR23126322.lite        11753888        2023-01-29 06:59:17     3554b2a57be48c519910a4faa4009f8c        1      SRA Lite Primary ETL     1       s3://sra-pub-zq-5/SRR23126322/SRR23126322.lite.1        s3.us-east-1    aws identity    https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR023/23126/SRR23126322/SRR23126322.lite.1    https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR023/23126/SRR23126322/SRR23126322.lite.1    worldwide       anonymous       gs://sra-pub-zq-102/SRR23126322/SRR23126322.lite.1     gs.us-east1      gcp identity    -       2021-12-03      cacao field     cacao tree      plant leaf      USA: Miami, FL  Theobroma cacao 25.64 N 80.29 W m1_2_1  MIGS/MIMS/MIMARKS.plant-associated      -       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/022/SRR23126322/SRR23126322_1.fastq.gz       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/022/SRR23126322/SRR23126322_2.fastq.gz       -       era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR231/022/SRR23126322/SRR23126322_1.fastq.gz    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR231/022/SRR23126322/SRR23126322_2.fastq.gz
# SRR23126323     SRP418302       Cacao phyllosphere microbiome is related to Phytophthora palmivora resistance   SRX19077291     amplicon sequencing of Theobroma cacao: phyllosphere microbiome amplicon sequencing of Theobroma cacao: phyllosphere microbiome 662107  phyllosphere metagenome m1_1_9_16S      AMPLICON        METAGENOMIC     PCR     PAIRED  SRS16492251     -       SAMN32798684    PRJNA925518     Illumina NovaSeq 6000  Illumina NovaSeq 6000    ILLUMINA        169266  17716421        169266  75999907        m1_1_9_1.fq     SRR23126323.lite        13123586        2023-01-29 06:59:16     c3ed4207699a6eb543d00e642f38dc35        1      SRA Lite Primary ETL     1       s3://sra-pub-zq-5/SRR23126323/SRR23126323.lite.1        s3.us-east-1    aws identity    https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR023/23126/SRR23126323/SRR23126323.lite.1    https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR023/23126/SRR23126323/SRR23126323.lite.1    worldwide       anonymous       gs://sra-pub-zq-102/SRR23126323/SRR23126323.lite.1     gs.us-east1      gcp identity    -       2021-11-19      cacao field     cacao tree      plant leaf      USA: Miami, FL  Theobroma cacao 25.64 N 80.29 W m1_1_9  MIGS/MIMS/MIMARKS.plant-associated      -       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/023/SRR23126323/SRR23126323_1.fastq.gz       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR231/023/SRR23126323/SRR23126323_2.fastq.gz       -       era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR231/023/SRR23126323/SRR23126323_1.fastq.gz    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR231/023/SRR23126323/SRR23126323_2.fastq.gz
```