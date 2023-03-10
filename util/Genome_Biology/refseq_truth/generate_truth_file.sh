# all commands used to generate the truth file

# create file to validate
USER_BINS_FILE="refseq.sample.ids"
READ_VALIDATION_FILE="read.validation.all"

grep -v QUERY_NAME ${USER_BINS_FILE} | cut -f 1 | awk 'BEGIN{OFS=":"}{print "dummy",$0}' > ${READ_VALIDATION_FILE}

# validate with yara
VALIDATION_RESULT="all_reads.validated"

./validate_all_reads.sh ${READ_VALIDATION_FILE} &> ${VALIDATION_RESULT}

# normalise result
READFILE="RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq";
VALIDATION_RESULT_SORTED="all_reads.validated.sorted";
QUERY_NAMES_FILE="query.names";

grep '@' ${READFILE} | tr --delete '@' | awk '{print $0":"}' | sort | tr --delete ':' > ${QUERY_NAMES_FILE}

# grep -v -w -f small.qnames all_reads.validated > all_reads.validated.250only
sort ${VALIDATION_RESULT} > ${VALIDATION_RESULT_SORTED}

NORMALISED_RESULT="the_one_and_only.normlised.truth"

# cpp file in util/applications/src/Genome_Biology
normalise_yara_truth_file --yara_results ${VALIDATION_RESULT_SORTED} --query_names ${QUERY_NAMES_FILE} --output_file ${NORMALISED_RESULT}

# create proper raptor file with user bin ids
RAPTOR_BASE_FILE="raptor.header"; # corresponds to the comparison.user_bin.ids file above
TRUTH_FILE="the_one_and_only.truth";

cat ${RAPTOR_BASE_FILE} > ${TRUTH_FILE}
cat ${NORMALISED_RESULT} >> ${TRUTH_FILE}
