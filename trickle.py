
from integrate import evolve
from initialise import initialise
from save_file import save_file
import matplotlib.pyplot as plt
from constants import DT_INIT, T_END, T_SAVE_EVERY
	
def main():

	data = initialise()

	dt = DT_INIT
	t = 0.
	t_end = T_END
	n_file = 0
	t_save_every = T_SAVE_EVERY
	step=0

	while t <= t_end:
		step += 1
		print("Step {}".format(step), "dt = {:.5E}".format(dt), "time = {:.5E}".format(t), "percentage = {:.2f}".format(t/T_END*100))

		data, dt = evolve(data, dt)

		t += dt
		  	
		if t >= n_file*t_save_every:
			
			print("Writing to disk file {}".format(n_file))
			save_file(data, n_file)

			plt.imshow(data[0], vmin=0.9, vmax=2.1)
			plt.savefig("Data/data0-{}.png".format(n_file))

			n_file+=1

if __name__ == "__main__":
	main()