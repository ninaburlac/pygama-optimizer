# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/lfs/l1/legend/users/dandrea/software/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/lfs/l1/legend/users/dandrea/software/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/lfs/l1/legend/users/dandrea/software/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/lfs/l1/legend/users/dandrea/software/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<


#python /lfs/l1/legend/users/dandrea/pygama/pgt/processing.py -r 117 -r2d -o > /lfs/l1/legend/users/dandrea/pygama/pgt/process_raw_to_dsp.out

python /lfs/l1/legend/users/dandrea/pygama/pygama-optimizer/optimizer.py -r 18 -c 2 -g -p -f -t -d /lfs/l1/legend/users/dandrea/pygama/pygama-optimizer > /lfs/l1/legend/users/dandrea/pygama/pygama-optimizer/optimizer-zac.out

