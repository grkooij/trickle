import h5py
import os

def save_file(data, n_file):

	filename = 'Data/'
	os.makedirs(os.path.dirname(filename), exist_ok=True)
		
	hf = h5py.File('Data/data.{}.hdf5'.format(n_file), 'w')

	hf.create_dataset('density', data=data[0])
	hf.create_dataset('vx', data=data[1])
	hf.create_dataset('vy', data=data[2])
	hf.create_dataset('pressure', data=data[3])

	hf.close()
	
	return