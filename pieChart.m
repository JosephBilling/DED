
clear
clc 
close all;

TPP=1;
p=1;
Floor=0;
Frames=20
a=1:500;
n=19


%Mp=linspace((n+1).*(TPP-n.*Floor)/(n.^2),0.19*(TPP-n.*Floor),Frames);
%Mp=linspace(0.1*(TPP-n.*Floor),0.20*(TPP-n.*Floor),Frames);
Mp=0.19*(TPP-n.*Floor)
j=1:n;
for k=1:length(Mp);     

    uNext=max(TPP,3);
    uNow=1;
    while(abs(uNow-uNext)>0.00000005)

        f=((TPP-n*Floor-Mp(k))*uNow^(n-p));
        df=(n-p)*((TPP-n*Floor-Mp(k))*uNow^(n-p-1));
        for i = 1:n
            if(i~=p)
                f=f-Mp(k)*uNow^(n-i);
                df=df-Mp(k)*(n-i)*uNow^(n-i-1);
            end
        end
        uNow=uNext
        uNext=uNow-f/df;
    end
    r(k)=log(uNext);
    A(k)=Mp(k).*exp((p-n).*r(k));
    R(k,j)=A(k).*exp((n-j).*r(k))+Floor;
    C(k,a)=A(k).*exp((n-a).*r(k))+Floor;
end
k


%pie3(R(length(Mp),j)

%pie3(ones([1 20]))
%title('Equal Resource Distribution')
figh = figure;

for k=1:length(Mp)
    
    clf
    
    Mp_k=Mp(k);
    R_k=R(k,j);
    %C_K=C(k,a);
    %A_K=A(k);
    %r_K=r(k);
   
        
    pie3(fliplr(R_k));
    %ylim([0 10]);
    %ylim([0 Mp(end)+Floor]);
    %xlim([0 n+1]);
    %hold on;
    %plot(a,C_K);
    %xlabel('Income Percentile','FontSize',8);
    %ylabel('Percent of total Income');
    title(['Exponential Resource Distribution when the tribe leader has ', num2str(100*Mp_k/TPP),' % of total resources'],'FontSize',8);
 
    
   % hold on;
   % a1=linspace(0,n+1,500)
   % Floora=Floor.*a1./a1
   % plot(a1,Floora);
  
    %legend('% Total Income',['y=',num2str(A_K),'*exp(-',num2str(r_K),'*x)'],['Floor',' = ',num2str(Floor),'%'])
    
    movieVector(k) = getframe(figh,[10 10 520 400]);
    
end

myWriter = VideoWriter('curve','MPEG-4');
myWriter.FrameRate = 70;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);
   