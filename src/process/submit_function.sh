#!/bin/bash
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true

submit_function() {
    local MAIN_DIR="$1"
    local FOLDER_NAME="$2"

    export SCRIPT="${MAIN_DIR}/src/process/main_process.sh"
    export WORK_DIR="${MAIN_DIR}/data/${FOLDER_NAME}"
    export LOG_DIR="${WORK_DIR}/logs"

    mkdir -p "${LOG_DIR}"
    echo -e "\n[INFO] Current working directory: $WORK_DIR"

    # For UGER
    if [[ $SERVER == "UGER" ]]; then
        TIME_HMS=$(printf "%02d:%02d:00" $((TIME/60)) $((TIME%60)))

        for sample in "${sample_list[@]}"; do

            out_log=${LOG_DIR}/${sample}.out
            out_err=${LOG_DIR}/${sample}.err
            rm -f $out_log
            rm -f $out_err

            if [[ $STEP == 0 ]]; then
                job_cmd="${SCRIPT} -w ${WORK_DIR} -ids ${sample} -s 1 -c ${CORE} -rp ${RERUN} && \
                            ${SCRIPT} -w ${WORK_DIR} -ids ${sample} -s 2 -c ${CORE} -rp ${RERUN} && \
                            ${SCRIPT} -w ${WORK_DIR} -ids ${sample} -s 3 -c ${CORE} -rp ${RERUN} && \
                            ${SCRIPT} -w ${WORK_DIR} -ids ${sample} -s 4 -c ${CORE} -rp ${RERUN} && \
                            ${SCRIPT} -w ${WORK_DIR} -ids ${sample} -s 5 -c ${CORE} -rp ${RERUN}"

                qsub -N "${sample}" -cwd \
                    -pe smp "$CORE" \
                    -l h_vmem="$MEM" \
                    -l h_rt="$TIME_HMS" \
                    -o "$out_log" \
                    -e "$out_err" \
                    -b y "$job_cmd"
            else
                qsub -N "${sample}" -cwd \
                    -pe smp "$CORE" \
                    -l h_vmem="$MEM" \
                    -l h_rt="$TIME_HMS" \
                    -o "$out_log" \
                    -e "$out_err" \
                    -b y "${SCRIPT} -w ${WORK_DIR} -ids ${sample} -s ${STEP} -c ${CORE} -rp ${RERUN}"
            fi
        done

    # For local
    elif [[ $SERVER == "LOCAL" ]]; then

        run_group() {
            local group_ids="$1"
            local group_name="$2"
            local log_file="${LOG_DIR}/log_batch_local_${group_name}.out"
            rm -f "$log_file"

            echo "Running: ${group_ids}"
            echo "Log: ${log_file}"

            export SCRIPT WORK_DIR CORE RERUN group_ids

            if [[ $STEP == 0 ]]; then
                nohup bash -c "
                    \"$SCRIPT\" -w \"$WORK_DIR\" -ids \"$group_ids\" -s 1 -c \"$CORE\" -rp \"$RERUN\" &&
                    \"$SCRIPT\" -w \"$WORK_DIR\" -ids \"$group_ids\" -s 2 -c \"$CORE\" -rp \"$RERUN\" &&
                    \"$SCRIPT\" -w \"$WORK_DIR\" -ids \"$group_ids\" -s 3 -c \"$CORE\" -rp \"$RERUN\" &&
                    \"$SCRIPT\" -w \"$WORK_DIR\" -ids \"$group_ids\" -s 4 -c \"$CORE\" -rp \"$RERUN\" &&
                    \"$SCRIPT\" -w \"$WORK_DIR\" -ids \"$group_ids\" -s 5 -c \"$CORE\" -rp \"$RERUN\"
                " > "$log_file" 2>&1 &
            else
                nohup bash -c "
                    \"$SCRIPT\" -w \"$WORK_DIR\" -ids \"$group_ids\" -s \"$STEP\" -c \"$CORE\" -rp \"$RERUN\"
                " > "$log_file" 2>&1 &
            fi
        }

        if [[ ${#sample_list[@]} -eq 1 ]]; then
            run_group "${sample_list[0]}" "${sample_list[0]}"
        else
            BATCH_SIZE=$(( ${#sample_list[@]} / 2 ))
            for ((i=0; i<${#sample_list[@]}; i+=BATCH_SIZE)); do
                group=("${sample_list[@]:i:BATCH_SIZE}")
                group_ids=$(IFS=, ; echo "${group[*]}")
                group_name=$(IFS=_ ; echo "${group[*]}")
                run_group "$group_ids" "$group_name"
            done
        fi

    fi

}