%% DED (Discrete Exponential Distribution) Tool
% The DED Tool takes inputs (Total Prize Pool) TPP, number of discrete points (n), and 
clear
clc 
close all;

%% Input Paramters
TPP=1376579; %Total Prize Pool - For every distrubtion, summing the value of x[n] for all n will equal TPP
n=10 %The number of dicrete points for the graph. 

%% Maximum value of x[1]/TPP. The DED tool creates a distrubtion for every value of Mp on its range and plots it as a function of time.
% The lower bound for Mp is TPP/n while the upper bound 
MpMax=0.50 



Frames=50 %Frames for animation

aMax=50; %Precision degree of the xplot
HappiHeight=0.15

error=0.00001 %Maximum allowable error for approximated exponential rate r(a)



Floor=1000


Mp=linspace((n+1)/(n^2)*(TPP-n.*Floor),MpMax*(TPP-n*Floor),Frames)
    %(n+1)/(n^2) ~ (1/n) such that Mp starts near its lower bound & (n+1)/(n^2) > (1/n) such that the approximation functions for n >> 1.
j=1:n
a=1:aMax

%x plot
%% Loop which calulates an exponential rate and base for each Frame
for Frame=1:Frames;     
  
    uNext=max(TPP,3);
    uNow=1;
    while(abs(uNow-uNext)>error)

        f=((TPP-n*Floor-Mp(Frame))*uNow^(n-1));
        df=(n-1)*((TPP-n*Floor-Mp(Frame))*uNow^(n-1-1));
        for i = 1:n
            if(i~=1)
                f=f-Mp(Frame)*uNow^(n-i);
                df=df-Mp(Frame)*(n-i)*uNow^(n-i-1);
            end
        end
        uNow=uNext;
        uNext=uNow-f/df;
    end
    r(Frame)=log(uNext);
    %Exponential rate for the curve on each frame
    B(Frame)=Mp(Frame).*exp((1-n).*r(Frame));
    %Exponential base for the curve on each frame
    x(Frame,j)=B(Frame).*exp((n-j).*r(Frame))+Floor;
    %Array of Discrete Points for each frame
    C(Frame,a)=B(Frame).*exp(n*(1-(a/aMax)).*r(Frame))+Floor;
    %Exponential curve for each frame
end

%%
figh = figure;
a1=linspace(0,n+1,500)
Floora=Floor.*a1./a1
Basea=TPP*MpMax*HappiHeight*0.223*log10(TPP/n).*a1./a1
for Frame=1:length(Mp)
    
    clf
    
    Mp_FRAME=Mp(Frame);
    x_FRAME=x(Frame,j);
    C_FRAME=C(Frame,a);
    B_FRAME=B(Frame);
    r_FRAME=r(Frame);
    
    happiness=0.223*log10(x_FRAME)
    stem(j,TPP*MpMax*HappiHeight*happiness)
    hold on;
    
    stem(j,x_FRAME,':diamondr','MarkerSize',3);
    ylim([0 10]);
    if(Floor<0)
        ylim([1.2*Floor Mp(end)+Floor]);
    else
        ylim([0 Mp(end)+Floor]);
    end
    xlim([0 n+1]);
    hold on;
    plot(a/max(a)*n,C_FRAME);
    xlabel('Income Percentile','FontSize',8);
    ylabel('Percent of total Income');
    title(['Exponential Income Distribution when the 1st percentile has ', num2str(100*Mp_FRAME/TPP),' % of total income'],'FontSize',8);
    %annotation('textbox', [0.4, 0.6, 0.1, 0.1], 'String',['Total Income = ',num2str(sum(R_k))])
    annotation('textbox', [0.4, 0.625, 0.1, 0.1], 'String',['Approximation error = ',num2str(100*((sum(x_FRAME)/TPP)-1)),'%'])
    
       
    
    annotation('textbox', [0.4, 0.55, 0.1, 0.1], 'String',['Average Happiness = ',num2str(sum(happiness)/n)])
    if(mod(n,2))
        annotation('textbox', [0.4, 0.475, 0.1, 0.1], 'String',['Median Happiness = ',num2str(happiness(n/2+0.5))])        
    else
        annotation('textbox', [0.4, 0.475, 0.1, 0.1], 'String',['Median Happiness = ',num2str(0.5*(happiness(n/2)+happiness(n/2+1)))])        
    end
    annotation('textbox', [0.4, 0.3, 0.1, 0.1], 'String',['Relative Happiness = 0.223*log10(Income)'])

    
    hold on;
    a1=linspace(0,n+1,500)
    plot(a1,Floora);
    hold on;
    plot(a1,Basea);
  
    legend('% Total Income',['y=',num2str(B_FRAME),'*exp(-',num2str(r_FRAME),'*x)'],['Floor',' = $',num2str(Floor)],['Equal Distribution Happiness = ',num2str(0.223*log10(TPP/n))])
    
    movieVector(Frame) = getframe(figh,[10 10 520 400]);
    drawnow
end

myWriter = VideoWriter('curve','MPEG-4');
myWriter.FrameRate = 25;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);
   