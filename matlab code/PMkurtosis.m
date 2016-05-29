clear
jk=4;
% Alldata 1 Time, 2 PM2.5 3 PM10
tic
load CityLocation.mat
Cityname=city.name;
CityLon=city.lon;
CityLat=city.lat;
BINRL=[1:0.1:4];
NB=length(BINRL);
Ku=zeros(4,NB);
NC=length(Cityname);


for i=1:NC-1
    C1=Cityname(i,:);
    xi=find(C1==' ');
    C1(xi)=[];

    load(['../source_data_mat/',C1,'AQI.mat']);
    T1=Alldata(:,1);
    PM1=Alldata(:,jk); %change index here
    NT=length(T1);
    if NT<100 % if the data length less than 100 skip
        continue;
    end

    for j=i+1:NC

        C2=Cityname(j,:);
        xi=find(C2==' ');
        C2(xi)=[];
        load(['../source_data_mat/',C2,'AQI.mat']);
        T2=Alldata(:,1);
        PM2=Alldata(:,jk); % change idnex here
        NT=length(T2);
        if NT<100
            continue;
        end

        dr=distance(CityLon(i),CityLat(i),CityLon(j),CityLat(j))*earthRadius('km')*pi/180;
        TS=max([T1(1) T2(1)]);
        TE=min([T1(end) T2(end)]);

        xi1=find(T1>=TS & T1<=TE);
        NX1=length(xi1);

        xi2=find(T2>=TS & T2<=TE);
        NX2=length(xi2);

        NX=min([NX1,NX2]);
        DPM=zeros(1,length(NX));

        TT=[T1' T2'];
        [T,I]=unique(TT);
        NT=length(T);
        k=0;
        for kk=1:NT
            xi1=find(T1==T(kk));
            xi2=find(T2==T(kk));

            if ~isempty(xi1) && ~isempty(xi2)
               k=k+1;
               DPM(k)=PM1(xi1(1))-PM2(xi2(1));
            end

        end
        DPM=(DPM(1:k));
%         DPM=DPM(DPM>0);

        dr=log10(dr);
        kk=fix((dr-1-0.05)/0.1+1);
        if kk>0 && kk<=NB
            Ku(1,kk)=Ku(1,kk)+nansum(DPM.^2);
            Ku(2,kk)=Ku(2,kk)+nansum(DPM.^3);
            Ku(3,kk)=Ku(3,kk)+nansum(DPM.^4);
            Ku(4,kk)=Ku(4,kk)+length(DPM(DPM>-1e5));
        end
    end
    toc
    i
end
Ku=bsxfun(@times,Ku,Ku(4,:).^-1);
save PM25kurtosis Ku BINRL
