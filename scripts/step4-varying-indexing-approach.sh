seed=0
result_dir=../data/result/step4

mkdir -p ${result_dir}

datasets=( aids yago )
methods=( alleyTPI )
ps=( 0.001 )

for data in ${datasets[@]}; do
    for method in ${methods[@]}; do
        export GCARE_ALLEY_TPI_NUM_GROUP=32

        # naive, maxL = 5
        export GCARE_ALLEY_TPI_FAIL_THRESHOLD=0.0
        export GCARE_ALLEY_TPI_MAXL=5

        ./build.sh ${method} ${data} ${seed} ${result_dir} 
        for p in ${ps[@]}; do
            ./est.sh ${method} ${data} ${p} ${seed} ${result_dir}
        done

        # naive, maxL = 4
        export GCARE_ALLEY_TPI_FAIL_THRESHOLD=0.0
        export GCARE_ALLEY_TPI_MAXL=4

        ./build.sh ${method} ${data} ${seed} ${result_dir} 
        for p in ${ps[@]}; do
            ./est.sh ${method} ${data} ${p} ${seed} ${result_dir}
        done

        # TPI, maxL = 5
        export GCARE_ALLEY_TPI_FAIL_THRESHOLD=0.9
        export GCARE_ALLEY_TPI_MAXL=5

        ./build.sh ${method} ${data} ${seed} ${result_dir} 
        for p in ${ps[@]}; do
            ./est.sh ${method} ${data} ${p} ${seed} ${result_dir}
        done

        # TPI, maxL = 4
        export GCARE_ALLEY_TPI_FAIL_THRESHOLD=0.9
        export GCARE_ALLEY_TPI_MAXL=4

        ./build.sh ${method} ${data} ${seed} ${result_dir} 
        for p in ${ps[@]}; do
            ./est.sh ${method} ${data} ${p} ${seed} ${result_dir}
        done
    done
done
