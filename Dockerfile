FROM broadinstitute/gatk:4.5.0.0

RUN apt update && apt install -y zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev

COPY install_samtools_tools.sh /tmp
RUN bash /tmp/install_samtools_tools.sh

# Download Picard
RUN mkdir /app && wget -P /app https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar

# Install BWA
RUN wget -P /app https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
RUN tar xjvf /app/bwa-0.7.17.tar.bz2 --directory /app
# fix for gcc 11.4.0 https://www.biostars.org/p/9521965/
RUN cd /app/bwa-0.7.17 && sed -i 's/^const uint8_t/extern const uint8_t/g' rle.h && make && ln bwa /bin

# Download references
RUN mkdir /references && wget -P /references https://storage.googleapis.com/genomics-public-data/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.fasta \
    https://storage.googleapis.com/genomics-public-data/references/hg38/v0/chrM/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta

# Prepare references
RUN bwa index /references/Homo_sapiens_assembly38.chrM.fasta \
    && bwa index /references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta \
RUN samtools faidx /references/Homo_sapiens_assembly38.chrM.fasta \
    samtools faidx /references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta
RUN java -jar /app/picard.jar CreateSequenceDictionary -R /references/Homo_sapiens_assembly38.chrM.fasta  \
    && java -jar /app/picard.jar CreateSequenceDictionary -R /references/Homo_sapiens_assembly38.chrM.shifted_by_8000_bases.fasta