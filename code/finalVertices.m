%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code written by Mauricelle and Vargas
% Last update: Jan 18, 2024
% Motivation: experimental data collected
% from a shaking table. Procedure that computes
% the vertices from a large database of matrices.
% E-mail: avargas@utfpr.edu.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all, close all, clc, format long, format compact,

disp(' .... procedure that computes the vertices from a large database (it may take some hours) ...')

fid = fopen('listaData.txt');
tline = fgetl(fid);
count = 1;
while ischar(tline)
    nome{count} = sprintf('%s',tline);
    tline = fgetl(fid);
    count = count+1;
end

fclose(fid);

A_aggreg = []; B_aggreg = [];
for cx=1:max(size(nome))
    text_file = sprintf('clean_matrices_%0.3i.mat',cx);
    load(text_file);
    A_aggreg = [A_aggreg  A_po];
    B_aggreg = [B_aggreg  B_po];
end

[A_vertices,B_vertices] = PolytopeGenerator(A_aggreg,B_aggreg);
savefile = sprintf('vertices_final.mat');
save(savefile, 'A_vertices', 'B_vertices','-v7');
