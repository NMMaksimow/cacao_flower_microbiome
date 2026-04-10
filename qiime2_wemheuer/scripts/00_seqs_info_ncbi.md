Поиск данных секвенирования ассоциированных с проектом PRJNA594470 через CLI tools
1. Метаданные через Entrez Direct

# Основная таблица всех ранов — размеры, имена, платформа, стратегия
```bash
esearch -db sra -query "PRJNA594470" | \
    efetch -format runinfo > PRJNA594470_runinfo.csv
```

# Посмотреть: Run, size_MB, LibraryStrategy, LibraryLayout, spots, SampleName и т.д.
```bash
column -t -s',' PRJNA594470_runinfo.csv | less -S
```

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

```bash
esearch -db sra -query "PRJNA594470" | efetch -format runinfo > runinfo.csv

head -3 runinfo.csv
# Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash
# SRR10614063,2020-03-12 09:40:51,2019-12-10 04:53:41,73219,43605277,73219,595,29,,https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614063.sralite.1,SRX7293343,Fun_T32,AMPLICON,PCR,METAGENOMIC,PAIRED,0,0,ILLUMINA,Illumina MiSeq,SRP235361,PRJNA594470,,594470,SRS5786637,SAMN13519336,simple,1297858,leaf metagenome,Tree32,,,,,,,no,,,,,GENOMIC AND APPLIED MICROBIOLOGY & GOETTINGEN GENOMICS LABORATORY,SRA1008811,,public,3678FF6290A9F27D1C3FBB95F3B0D8EF,0EC9187FFA58135EF4E869DBE14AAB35
# SRR10614062,2020-03-12 09:40:51,2019-12-10 04:53:16,123636,74168797,123636,599,46,,https://sra-downloadb.be-md. ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614062.sralite.1,SRX7293344,Fun_T33,AMPLICON,PCR,METAGENOMIC,PAIRED,0,0,ILLUMINA,Illumina MiSeq,SRP235361,PRJNA594470,,594470,SRS5786638,SAMN13519337,simple,1297858,leaf metagenome,Tree33,,,,,,,no,,,,,GENOMIC AND APPLIED MICROBIOLOGY & GOETTINGEN GENOMICS LABORATORY,SRA1008811,,public,CDA6478EC7406E77FA31F5169A225B40,C32BD2E8DEA4869D9F76BD5AE690AE14

wc -l runinfo.csv 
# 129 runinfo.csv

awk -F',' 'NR>1 {sum+=$8} END {printf "Total: %.1f GB\n", sum/1024}' runinfo.csv
# Total: 3.7 GB

pysradb metadata --detailed PRJNA594470 > metadata.tsv
head -3 metadata.tsv 
# run_accession   study_accession study_title     experiment_accession    experiment_title        experiment_desc organism_taxid  organism_name   library_name    library_strategy        library_source  library_selection       library_layoutsample_accession        sample_title    biosample       bioproject      instrument      instrument_model        instrument_model_desc   total_spots     total_size      run_total_spots run_total_bases run_alias       public_filename public_size   public_date     public_md5      public_version  public_semantic_name    public_supertype        public_sratoolkit       aws_url aws_free_egress aws_access_type public_url      ncbi_url        ncbi_free_egress        ncbi_access_type      gcp_url gcp_free_egress gcp_access_type experiment_alias        collection_date env_broad_scale env_local_scale env_medium      geo_loc_name    host    lat_lon altitude        temp    tree_height     plantation      replicate    biosamplemodel   ena_fastq_http  ena_fastq_http_1        ena_fastq_http_2        ena_fastq_ftp   ena_fastq_ftp_1 ena_fastq_ftp_2
# SRR10614059     SRP235361       Theobroma cacao leaf metagenome SRX7293347      Fungal endophyte communities in leaves of Theobroma cacao, sample Tree36        Fungal endophyte communities in leaves of Theobroma cacao, sample Tree36     1297858  leaf metagenome Fun_T36 AMPLICON        METAGENOMIC     PCR     PAIRED  SRS5786641      -       SAMN13519340    PRJNA594470     Illumina MiSeq  Illumina MiSeq  ILLUMINA        76746   29965293        76746   45406069        Fun_T36_R1.fastq.gz   SRR10614059.sralite     11694441        2020-06-24 01:50:39     798c2a53ffd3a5a32dcabcf015ce2786        1       SRA Lite        Primary ETL     1       s3://sra-pub-zq-8/SRR10614059/SRR10614059.sralite.1     s3.us-east-1 aws identity     https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR010/614/SRR10614059.sralite.1        https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos5/sra-pub-zq-14/SRR010/614/SRR10614059.sralite.1        worldwide    anonymous        gs://sra-pub-zq-108/SRR10614059/SRR10614059.zq.1        gs.us-east1     gcp identity    -       2014    cropland biome  cacao plantation        leaf    Cameroon: Bakoa Theobroma cacao 4.5739 N 11.1798 E      683m    24.878.3m     Bakoa   Bakoa_rep4      MIGS/MIMS/MIMARKS.plant-associated      -       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/059/SRR10614059/SRR10614059_1.fastq.gz       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/059/SRR10614059/SRR10614059_2.fastq.gz     -       era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/059/SRR10614059/SRR10614059_1.fastq.gz    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/059/SRR10614059/SRR10614059_2.fastq.gz
# SRR10614060     SRP235361       Theobroma cacao leaf metagenome SRX7293346      Fungal endophyte communities in leaves of Theobroma cacao, sample Tree35        Fungal endophyte communities in leaves of Theobroma cacao, sample Tree35     1297858  leaf metagenome Fun_T35 AMPLICON        METAGENOMIC     PCR     PAIRED  SRS5786640      -       SAMN13519339    PRJNA594470     Illumina MiSeq  Illumina MiSeq  ILLUMINA        10164   4085271 10164   5967227 Fun_T35_R1.fastq.gz  SRR10614060.sralite      1606168 2020-06-20 10:57:53     a9fb63a3500fe50a5e6fa259d7b54c86        1       SRA Lite        Primary ETL     1       s3://sra-pub-zq-7/SRR10614060/SRR10614060.sralite.1     s3.us-east-1    aws identity    https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614060.sralite.1     https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos9/sra-pub-zq-922/SRR010/614/SRR10614060.sralite.1       worldwide       anonymous       gs://sra-pub-zq-108/SRR10614060/SRR10614060.zq.1      gs.us-east1     gcp identity    -       2014    cropland biome  cacao plantation        leaf    Cameroon: Bakoa Theobroma cacao 4.5739 N 11.1798 E      683m    24.87   7.8m    Bakoa   Bakoa_rep3    MIGS/MIMS/MIMARKS.plant-associated      -       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/060/SRR10614060/SRR10614060_1.fastq.gz       http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/060/SRR10614060/SRR10614060_2.fastq.gz       -    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/060/SRR10614060/SRR10614060_1.fastq.gz     era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/SRR106/060/SRR10614060/SRR10614060_2.fastq.gz
```

