import dpdata
import numpy as np

vasp_multi_systems = dpdata.MultiSystems.from_dir(
    dir_name="./", file_name="*OUTCAR", fmt="vasp/outcar"
)
training_set = dpdata.MultiSystems()
validation_set = dpdata.MultiSystems()
#print(vasp_multi_systems.systems)
#vasp_multi_systems.to_deepmd_raw("./deepmd_data/")
#vasp_multi_systems.to_deepmd_npy("./deepmd_data/")

for system in vasp_multi_systems:
    frame_num = len(system)
    index_validation = np.random.choice(frame_num, size=int(0.2*frame_num), replace=False)
    index_training = list(set(range(frame_num))-set(index_validation))

    data_training = system.sub_system(index_training)
    training_set.append(data_training)
    if len(index_validation)>0:
        data_validation = system.sub_system(index_validation)
        validation_set.append(data_validation)
#    print('# the training data contains %d frames' % len(index_training))
#    print('# the validation data contains %d frames' % len(index_validation))

print(training_set)
print(validation_set)
training_set.to_deepmd_raw('training_data')
training_set.to_deepmd_npy('training_data')
validation_set.to_deepmd_raw('validation_data')
validation_set.to_deepmd_npy('validation_data')
