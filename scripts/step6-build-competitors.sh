seed=0
result_dir=../data/result/step6

mkdir -p ${result_dir}

datasets=( lubm80 aids human yago hprd youtube )
methods=( cset sumrdf bsk )

for data in ${datasets[@]}; do
    for method in ${methods[@]}; do
        
        # set method-specific parameters
        if [ ${method} == bsk ] # Bound sketch requires one more parameter: the budget of BSK
        then
            export GCARE_BSK_BUDGET=4096
        fi

        if [ ${method} == sumrdf ]  # SumRDF requires one more parameter: the threshold to determine whether merging two summary vertices
        then
            export GCARE_SUMRDF_THRESHOLD=0.15
        fi

        if [ ${method} == alleyTPI ]; then
            export GCARE_ALLEY_TPI_FAIL_THRESHOLD=0.9
            export GCARE_ALLEY_TPI_NUM_GROUP=32
            if [ ${data} == hprd ] || [ ${data} == human ] || [ ${data} == youtube ]; then
                export GCARE_ALLEY_TPI_MAXL=4
            else
                export GCARE_ALLEY_TPI_MAXL=5
            fi
        fi

				./build.sh ${method} ${data} ${seed} ${result_dir} 
    done
done
