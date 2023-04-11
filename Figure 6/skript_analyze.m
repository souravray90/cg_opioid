%% Author: mwema <mwema@LAPTOP-BNO5OPFS>
%% Created: 2022-07-16
clear;

ph7_3f_prot=dlmread('com_prot_7_3f.xvg', '', 26,0);
ph7_3f_prot=[ph7_3f_prot;ph7_3f_prot;ph7_3f_prot;ph7_3f_prot];
ph7_3f=dlmread('com_lig1_7_3f.xvg', '', 26,0);
ph7_3f=[ph7_3f;dlmread('com_lig2_7_3f.xvg', '', 26,0)];
ph7_3f=[ph7_3f;dlmread('com_lig3_7_3f.xvg', '', 26,0)];
ph7_3f=[ph7_3f;dlmread('com_lig4_7_3f.xvg', '', 26,0)];

ph7_f_prot=dlmread('com_prot_7_f_p.xvg', '', 26,0);
ph7_f_prot=[ph7_f_prot;ph7_f_prot;ph7_f_prot;ph7_f_prot];
ph7_f=dlmread('com_lig1_7_f_p.xvg', '', 26,0);
ph7_f=[ph7_f;dlmread('com_lig2_7_f_p.xvg', '', 26,0)];
ph7_f=[ph7_f;dlmread('com_lig3_7_f_p.xvg', '', 26,0)];
ph7_f=[ph7_f;dlmread('com_lig4_7_f_p.xvg', '', 26,0)];

ph5_3f_prot=dlmread('com_prot_5_3f_p_dsb.xvg', '', 26,0);
ph5_3f_prot=[ph5_3f_prot;ph5_3f_prot;ph5_3f_prot;ph5_3f_prot];
ph5_3f=dlmread('com_lig1_5_3f_p_dsb.xvg', '', 26,0);
ph5_3f=[ph5_3f;dlmread('com_lig2_5_3f_p_dsb.xvg', '', 26,0)];
ph5_3f=[ph5_3f;dlmread('com_lig3_5_3f_p_dsb.xvg', '', 26,0)];
ph5_3f=[ph5_3f;dlmread('com_lig4_5_3f_p_dsb.xvg', '', 26,0)];

ph5_f_prot=dlmread('com_prot_5_f_p_dsb.xvg', '', 26,0);
ph5_f_prot=[ph5_f_prot;ph5_f_prot;ph5_f_prot;ph5_f_prot];
ph5_f=dlmread('com_lig1_5_f_p_dsb.xvg', '', 26,0);
ph5_f=[ph5_f;dlmread('com_lig2_5_f_p_dsb.xvg', '', 26,0)];
ph5_f=[ph5_f;dlmread('com_lig3_5_f_p_dsb.xvg', '', 26,0)];
ph5_f=[ph5_f;dlmread('com_lig4_5_f_p_dsb.xvg', '', 26,0)];

ph7_3f=ph7_3f-ph7_3f_prot;
%ph7_3f=ph7_3f(1:1001, :);
ph7_f=ph7_f-ph7_f_prot;
ph5_3f=ph5_3f-ph5_3f_prot;
ph5_f=ph5_f-ph5_f_prot;

figure(1); 
plot(ph7_f(:,2), ph7_f(:,3),'r.'); hold on;
plot(ph7_3f(:,2), ph7_3f(:,3),'b.'); hold off;
figure(7);
plot(ph5_f(:,2), ph5_f(:,3),'r.'); hold on;
plot(ph5_3f(:,2), ph5_3f(:,3),'b.'); hold off;

voronoi_centers=[];
for i=-4.5:0.5:4.5
  for j=-4.5:0.5:4.5
    voronoi_centers=[voronoi_centers;[i, j]];
  end
end
figure(3);
plot(voronoi_centers(:,1),voronoi_centers(:,2),'*'); 
  
figure(2);
[a,b]=hist(ph5_f(:,4),60 ); plot(b,a); hold on;
[a,b]=hist(ph7_f(:,4), 60); plot(b,a); 
[a,b]=hist(ph5_3f(:,4), 60); plot(b,a); 
[a,b]=hist(ph7_3f(:,4), 60); plot(b,a);  
 


count_7_f=zeros(length(voronoi_centers),1);
for i=1:4004
  dist=voronoi_centers-ph7_f(i,[2 3]);
  dist=dist.^2; dist=sum(dist,2);
  [val,index]=min(dist);
  count_7_f(index)=count_7_f(index)+1;
end

count_5_3f=zeros(length(voronoi_centers),1);
for i=1:4004
  dist=voronoi_centers-ph5_3f(i,[2 3]);
  dist=dist.^2; dist=sum(dist,2);
  [val,index]=min(dist);
  count_5_3f(index)=count_5_3f(index)+1;
end

count_5_f=zeros(length(voronoi_centers),1);
for i=1:4004
  dist=voronoi_centers-ph5_f(i,[2 3]);
  dist=dist.^2; dist=sum(dist,2);
  [val,index]=min(dist);
  count_5_f(index)=count_5_f(index)+1;
end

count_7_3f=zeros(length(voronoi_centers),1);
for i=1:1001 %4004
  dist=voronoi_centers-ph7_3f(i,[2 3]);
  dist=dist.^2; dist=sum(dist,2);
  [val,index]=min(dist);
  count_7_3f(index)=count_7_3f(index)+1;
end

admat=zeros(length(voronoi_centers), length(voronoi_centers));
for i=1:length(voronoi_centers)
  dist=voronoi_centers-voronoi_centers(i,:);
  dist=dist.^2;
  dist=sqrt(sum(dist,2));
  index=find(dist<0.55);
  admat(i,index)=1;
end
admat=admat-eye(length(voronoi_centers));

ind_7_3f=find(count_7_3f>0);
Q7_3f=diag(1./sqrt(count_7_3f(ind_7_3f)))*admat(ind_7_3f, ind_7_3f)*diag(sqrt(count_7_3f(ind_7_3f)));
Q7_3f=Q7_3f-diag(sum(Q7_3f,2));
ind_7=find(count_7_f>0);
Q7_f=diag(1./sqrt(count_7_f(ind_7)))*admat(ind_7,ind_7)*diag(sqrt(count_7_f(ind_7)));
Q7_f=Q7_f-diag(sum(Q7_f,2));

ind_5_3f=find(count_5_3f>0);
Q5_3f=diag(1./sqrt(count_5_3f(ind_5_3f)))*admat(ind_5_3f, ind_5_3f)*diag(sqrt(count_5_3f(ind_5_3f)));
Q5_3f=Q5_3f-diag(sum(Q5_3f,2));
ind_5=find(count_5_f>0);
Q5_f=diag(1./sqrt(count_5_f(ind_5)))*admat(ind_5, ind_5)*diag(sqrt(count_5_f(ind_5)));
Q5_f=Q5_f-diag(sum(Q5_f,2));

figure(4); hold on;
[u73,v73]=eig(Q7_3f);
[val, ind]=sort(-diag(v73)); 
val73=val(2);
[val, indmax]=max(count_7_3f(ind_7_3f));
u73=u73(:,ind); u73=u73*sign(u73(indmax,2)); u73=(max(u73(:,2))-u73(:,2))/(max(u73(:,2))-min(u73(:,2)));
plot3(voronoi_centers(ind_7_3f,1), voronoi_centers(ind_7_3f,2), u73,'*'); 

[u53,v53]=eig(Q5_3f);
[val, ind]=sort(-diag(v53)); 
val53=val(2);
[val, indmax]=max(count_5_3f(ind_5_3f));
u53=u53(:,ind); u53=u53*sign(u53(indmax,2)); u53=(max(u53(:,2))-u53(:,2))/(max(u53(:,2))-min(u53(:,2)));
plot3(voronoi_centers(ind_5_3f,1), voronoi_centers(ind_5_3f,2), u53,'*'); 
hold off;

figure(5); hold on;
[u7,v7]=eig(Q7_f);
[val, ind]=sort(-diag(v7)); 
val7=val(2);
[val, indmax]=max(count_7_f(ind_7));
u7=u7(:,ind); u7=u7*sign(u7(indmax,2)); u7=(max(u7(:,2))-u7(:,2))/(max(u7(:,2))-min(u7(:,2)));
plot3(voronoi_centers(ind_7,1), voronoi_centers(ind_7,2), u7,'*'); 

[u5,v5]=eig(Q5_f);
[val, ind]=sort(-diag(v5)); 
val5=val(2);
[val, indmax]=max(count_5_f(ind_5));
u5=u5(:,ind); u5=u5*sign(u5(indmax,2)); u5=(max(u5(:,2))-u5(:,2))/(max(u5(:,2))-min(u5(:,2)));
plot3(voronoi_centers(ind_5,1), voronoi_centers(ind_5,2), u5,'*'); 
hold off;

% Ligand-Ligand Distance
ph7_minlig_3f=[];
for i=1:1001
    mindist=100000;
    for j=1:3
        for k=j+1:4
            dist=norm(ph7_3f(i+(j-1)*1001,2:4)-ph7_3f(i+(k-1)*1001,2:4));
            if (dist < mindist)
                mindist=dist;
            end
        end
    end
    ph7_minlig_3f=[ph7_minlig_3f;mindist];
end
ph7_minlig=[];
for i=1:1001
    mindist=100000;
    for j=1:3
        for k=j+1:4
            dist=norm(ph7_f(i+(j-1)*1001,2:4)-ph7_f(i+(k-1)*1001,2:4));
            if (dist < mindist)
                mindist=dist;
            end
        end
    end
    ph7_minlig=[ph7_minlig;mindist];
end

ph5_minlig_3f=[];
for i=1:1001
    mindist=100000;
    for j=1:3
        for k=j+1:4
            dist=norm(ph5_3f(i+(j-1)*1001,2:4)-ph5_3f(i+(k-1)*1001,2:4));
            if (dist < mindist)
                mindist=dist;
            end
        end
    end
    ph5_minlig_3f=[ph5_minlig_3f;mindist];
end
ph5_minlig=[];
for i=1:1001
    mindist=100000;
    for j=1:3
        for k=j+1:4
            dist=norm(ph5_f(i+(j-1)*1001,2:4)-ph5_f(i+(k-1)*1001,2:4));
            if (dist < mindist)
                mindist=dist;
            end
        end
    end
    ph5_minlig=[ph5_minlig;mindist];
end

xy_dist_pH7_3f=[[1:4004]',sqrt(sum(ph7_3f(:,2:3).^2,2))];

xy_dist_pH7_3f_lig1=xy_dist_pH7_3f(1:1001,:);
xy_dist_pH7_3f_lig2=xy_dist_pH7_3f(1002:2002,:);
xy_dist_pH7_3f_lig3=xy_dist_pH7_3f(2003:3003,:);
xy_dist_pH7_3f_lig4=xy_dist_pH7_3f(3004:4004,:);

z_pH7_3f=[[1:4004]', ph7_3f(:,4)];
z_pH7_3f_lig1=z_pH7_3f(1:1001,:);
z_pH7_3f_lig2=z_pH7_3f(1002:2002,:);
z_pH7_3f_lig3=z_pH7_3f(2003:3003,:);
z_pH7_3f_lig4=z_pH7_3f(3004:4004,:);


save session_marcus.mat
