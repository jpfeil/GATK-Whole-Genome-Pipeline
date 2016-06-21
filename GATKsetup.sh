# -e for exit, -x for echo, -o pipefail for report failure if ANY process fails (default is to report success if first)
set -e
set -x 
set -o pipefail

#get reference genome
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.gz
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta.fai.gz
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.dict.gz
#gunzip ucsc.hg19.fasta.gz
#gunzip ucsc.hg19.fasta.fai.gz 
#gunzip ucsc.hg19.dict.gz

#phase2 reference
#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
#gunzip hs37d5.fa.gz
samtools faidx hs37d5.fa

# VCF files for VariantRecalibrator_INDEL.RealignerTargetCreator knowns, and dbsnp for BaseRecalibraton
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
gunzip 1000G_phase1.indels.hg19.sites.vcf.gz
gunzip dbsnp_138.hg19.vcf.gz

#Resource files for VariantRecalibrator_SNP
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
gunzip hapmap_3.3.hg19.sites.vcf.gz
gunzip 1000G_omni2.5.hg19.sites.vcf.gz
gunzip 1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz


cd ~
#TOOLS

#Java
sudo add-apt-repository -y ppa:webupd8team/java
echo debconf shared/accepted-oracle-license-v1-1 select true | sudo debconf-set-selections
echo debconf shared/accepted-oracle-license-v1-1 seen true | sudo debconf-set-selections
sudo apt-get update
sudo apt-get install oracle-java7-set-default

# install s3cmd
sudo apt-get -qy install s3cmd

# get/install samtools
sudo apt-get -qy install zlib1g-dev
sudo apt-get -qy install samtools

# get/install picard, requires unzip program
wget https://github.com/broadinstitute/picard/releases/download/1.130/picard-tools-1.130.zip
sudo apt-get install unzip
unzip picard-tools-1.130.zip

# get GATK.jar
wget https://s3-us-west-2.amazonaws.com/bd2k-artifacts/10k-exomes/GenomeAnalysisTK.jar

#cat > ~/.s3cfg <<-END
#	[default]
#	access_key = $ACCESS_KEY
#	bucket_location = US
#	cloudfront_host = cloudfront.amazonaws.com
#	default_mime_type = binary/octet-stream
#	delete_removed = False
#	dry_run = False
#	enable_multipart = True
#	encoding = UTF-8
#	encrypt = False
#	follow_symlinks = False
#	force = False
#	get_continue = False
#	gpg_command = /usr/bin/gpg
#	gpg_decrypt = %(gpg_command)s -d --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
#	gpg_encrypt = %(gpg_command)s -c --verbose --no-use-agent --batch --yes --passphrase-fd %(passphrase_fd)s -o %(output_file)s %(input_file)s
#	gpg_passphrase =
#	guess_mime_type = True
#	host_base = s3.amazonaws.com
#	host_bucket = %(bucket)s.s3.amazonaws.com
#	human_readable_sizes = False
#	invalidate_on_cf = False
#	list_md5 = False
#	log_target_prefix =
#	mime_type =
#	multipart_chunk_size_mb = 15
#	preserve_attrs = True
#	progress_meter = True
#	proxy_host =
#	proxy_port = 0
#	recursive = False
#	recv_chunk = 4096
#	reduced_redundancy = False
#	secret_key = $SECRET_KEY
#	send_chunk = 4096
#	simpledb_host = sdb.amazonaws.com
#	skip_existing = False
#	socket_timeout = 300
#	urlencoding_mode = normal
#	use_https = False
#	verbosity = WARNING
#	website_endpoint = http://%(bucket)s.s3-website-%(location)s.amazonaws.com/
#	website_error =
#	website_index = index.html
END
