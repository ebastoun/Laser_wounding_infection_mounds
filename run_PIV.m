% Script to run urapiv for successive frames
for k=2:107
urapiv('/Users/effiebastounis/Documents/MDCKs/My_paper/final_revision_2020_11_21/Cell_star_protocols_2021_01_26/scripts/Pos2/ecad/','/Users/effiebastounis/Documents/MDCKs/My_paper/final_revision_2020_11_21/Cell_star_protocols_2021_01_26/scripts/Pos2/urapiv/','ecad',1,32,16,1,2,0,1,0,0,0*[1 1 1 1],0,2,3,3,8,[k-1],[k],1);
end