%% Run repositionNeuronLayer on all cells in population

cell_ids=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
cell_layers=[1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5];
layer_set_num = 10; 
nrn_pop_name = 'nrn_pop1';
nrn_model_ver = 'maxH'; 
opts = struct();
opts.par_on = 0;  % set to 1 to parallelize over available CPUs
opts.mode = 'both'; % see repositionNeuronLayer.m for definition
for i = 1:length(cell_ids)
    repositionNeuronLayer(cell_ids(i),cell_layers(i),layer_set_num,...
                          nrn_pop_name,nrn_model_ver,opts)
end
%% Compile repositioned neuron population data and generate new NeuronPop object

compileRepositionedPop(layer_set_num,nrn_pop_name,nrn_model_ver,'mode',...
                        opts.mode)