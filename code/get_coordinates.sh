wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP322526&result=read_run&fields=sample_accession,location,isolate,sample_alias" > PRJNA731936_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP305791&result=read_run&fields=sample_accession,location,isolate,sample_alias" > PRJNA700635_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP272267&result=read_run&fields=sample_accession,location,isolate,sample_alias" > PRJNA646656_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP254726&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA616181_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP149293&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA473595_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP111200&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA393325_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP254726&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA616181_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP295118&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680553_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294059&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680573_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294011&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680571_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294016&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680566_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294033&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680562_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294041&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680554_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP294009&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680541_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP293991&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA680508_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP121945&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJEB38505_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP063647&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA275184_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP045078&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA256477_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP045075&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA256249_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP045073&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA256476_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP045059&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA256246_locations.txt
wget -q  -O - "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP361743&result=read_run&fields=sample_accession,location,isolate,sample_alias" > 	PRJNA795879_locations.txt

cat *_locations.txt > locations.txt