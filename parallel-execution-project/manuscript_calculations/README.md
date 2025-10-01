# TO RUN CALCULATIONS SEQUENTIALLY AS A BATCH (with multiprocessing turned on):

for suffix in 02 04 06 08 10; do
  python water_0${suffix}_qcmanybody.py &> water_0${suffix}.py || break
done

# TO CONVERT HMBE input files to QCManyBody files

python hmbe_to_qcmanybody.py water_006.inp --bsse-typ
e "nocp" --max-workers 6

# NOTE: See script for full set of flags