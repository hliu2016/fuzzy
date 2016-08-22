***README for hybrid modeling of signaling pathway ************
1.	Synthetic data
The synthetic time-series data was produced by using BoolNet on the toy pathway. 
Just run toy.m to reproduce the experiment on synthetic data. Toy.m will load the data saved in ts_synthetic.mat and call the objfun_parallel_ts.m for Genetic Algorithm. 
2.	DREAM4 dataset
Run dream4.m to execute the experiment on the DREAM4 signaling pathway dataset.
3.	Cell fate prediction
The cell fate data was published by the Yaffe¡¯s search group. Please run Yaffe.m to learn the signaling pathway based on the time-series phosphoproteomics data.
