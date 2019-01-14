function logerror(M,L)

n=L-1;
vec_arnold=cell(n,1);
vec_gs=cell(n,1);
vec_gs_disp=cell(n,1);
vec_gs_sigma=cell(n,1);

string1_arndol='Arnold_';
string1_gs='BlockGS_';
string1_gs_disp='BlockGS_disp';
string1_gs_sigma='BlockGS_sigma';
string3='.txt';
string2=num2str(M);
string3='_';
string5='.txt';

figure
hold all
legend('-DynamicLegend');
for lev=2:n
string4=num2str(lev);
filename=strcat(string1_arndol,string2,string3,string4,string5);
vec_arnold{lev-1}=load(filename);
filename=strcat(string1_gs,string2,string3,string4,string5);    
vec_gs{lev-1}=load(filename);
filename=strcat(string1_gs_sigma,string2,string3,string4,string5);    
vec_gs_sigma{lev-1}=load(filename);
filename=strcat(string1_gs_disp,string2,string3,string4,string5);    
vec_gs_disp{lev-1}=load(filename);
end

for lev=2:n+1
string4=num2str(lev);
filename=strcat(string1_arndol,string2,string3,string4,string5);
end

mark=['^','s','o','.','x'];
for lev=2:n+1
    plot(vec_arnold{lev-1},'r','Marker',mark(lev-1));
    string1='Arnold_';
    string2=num2str(lev);
    string{lev-1}=strcat(string1,string2);
end
legend(string)
figure
hold all
for lev=2:n
    plot(vec_gs{lev-1},'r','Marker',mark(lev-1));
    string1='BLockGS_';
    string2=num2str(lev);
    string{lev-1}=strcat(string1,string2);
end
legend(string)

figure
hold all
for lev=2:n
    plot(vec_gs_sigma{lev-1},'r','Marker',mark(lev-1));
    plot(vec_gs_disp{lev-1},'g','Marker',mark(lev-1));
    string1_sigma='BLockGS_{sigma}L=';
    string1_disp='BLockGS_{disp}L=';
    string2=num2str(lev);
    string{2*(lev-1)-1}=strcat(string1_sigma,string2);
    string{2*(lev-1)}=strcat(string1_disp,string2);
end
legend(string)
% plot(vec_gs{lev-1});
% plot(vec_gs_sigma{lev-1});
% plot(vec_gs_disp{lev-1});
end