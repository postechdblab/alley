method=$1
data=$2
p=$3
seed=$4
result_dir=$5

input_dir=../data/queryset/${data}/
data_dir=../data/dataset/${data}/

if [ ${method} == cset ] || [ ${method} == bsk ] || [ ${method} == sumrdf ] # deterministic estimator
then
    repeat=1
else
    repeat=30
fi


# choose binary
if [ ${method} == cs ] || [ ${method} == bsk ]
then
    bin=../build/gcare_relation
elif [ ${method} == alley ] || [ ${method} == alleyTPI ]
then
    bin=../build/gcare_inter
else
    bin=../build/gcare_graph
fi


# set output location
if [ ${method} == bsk ]
then
    output_file=${result_dir}/${data}_${method}_b${GCARE_BSK_BUDGET}_s${seed}_query_result.txt
elif [ ${method} == sumrdf ]
then
    output_file=${result_dir}/${data}_${method}_p${p}_t${GCARE_SUMRDF_THRESHOLD}_s${seed}_query_result.txt
elif [ ${method} == alleyTPI ]
then
    output_file=${result_dir}/${data}_${method}_p${p}_L${GCARE_ALLEY_TPI_MAXL}_g${GCARE_ALLEY_TPI_NUM_GROUP}_t${GCARE_ALLEY_TPI_FAIL_THRESHOLD}_s${seed}_query_result.txt
else
    output_file=${result_dir}/${data}_${method}_p${p}_s${seed}_query_result.txt
fi


cmd="${bin} -q -m ${method} -i ${input_dir} -d ${data_dir}/${data} -p ${p} -n ${repeat} -s ${seed} -o ${output_file}"
${cmd} > "${output_file%.*}".log
