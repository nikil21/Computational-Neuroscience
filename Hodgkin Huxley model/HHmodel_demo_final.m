%THIS PROGRAM DEMONSTRATES HODGKIN HUXLEY MODEL IN CURRENT CLAMP EXPERIMENTS AND SHOWS ACTION POTENTIAL PROPAGATION
%Time is in secs, voltage in mvs, conductances in m mho/mm^2, capacitance in uF/mm^2

% threshold value of current is 0.0223
k=1;
i_step_size=0.0001;
for ImpCur=0:i_step_size:0.6
%ImpCur=input('enter the value of the impulse current in microamperes: ');
%TimeTot=input('enter the time for which stimulus is applied in milliseconds');

gkmax=.36;
vk=-77; 
gnamax=1.20;
vna=50; 
gl=0.003;
vl=-54.387; 
cm=.01; 

dt=0.01;
niter=100000;
t=(1:niter)*dt;
iapp=ImpCur*ones(1,niter);
%for i=1:100
 %   iapp(1,i)=ImpCur;
 %end;
v=-64.9964;
m=0.0530;
h=0.5960;
n=0.3177;

gnahist=zeros(1,niter);
gkhist=zeros(1,niter);
vhist=zeros(1,niter);
mhist=zeros(1,niter);
hhist=zeros(1,niter);
nhist=zeros(1,niter);


for iter=1:niter
  gna=gnamax*m^3*h; 
  gk=gkmax*n^4; 
  gtot=gna+gk+gl;
  vinf = ((gna*vna+gk*vk+gl*vl)+ iapp(iter))/gtot;
  tauv = cm/gtot;
  v=vinf+(v-vinf)*exp(-dt/tauv);
  alpham = 0.1*(v+40)/(1-exp(-(v+40)/10));
  betam = 4*exp(-0.0556*(v+65));
  alphan = 0.01*(v+55)/(1-exp(-(v+55)/10));
  betan = 0.125*exp(-(v+65)/80);
  alphah = 0.07*exp(-0.05*(v+65));
  betah = 1/(1+exp(-0.1*(v+35)));
  taum = 1/(alpham+betam);
  tauh = 1/(alphah+betah);
  taun = 1/(alphan+betan);
  minf = alpham*taum;
  hinf = alphah*tauh;
  ninf = alphan*taun;
  m=minf+(m-minf)*exp(-dt/taum);
  h=hinf+(h-hinf)*exp(-dt/tauh);
  n=ninf+(n-ninf)*exp(-dt/taun);
  vhist(iter)=v; mhist(iter)=m; hhist(iter)=h; nhist(iter)=n;
end
j=1;
realpeaks=zeros;
[peaks, locs]=findpeaks(vhist);
for temp=1:length(peaks) 
    if   peaks(temp) >=10 % minimum value at which a waveform is considered AP.
        realpeaks(j)=peaks(temp);        
        j=j+1;
    end
end
if realpeaks ~= 0
no_of_peaks(k)=length(realpeaks);
else
    no_of_peaks(k)=0;
end
k=k+1;
end
figure(1)
%subplot(2,1,1)
plot(t,vhist)
title('voltage vs time')

figure(2)
%subplot(2,1,2)
plot(t,mhist,'y-', t,hhist,'g.',t,nhist,'b-')
legend('m','h','n')

figure(3)
gna=gnamax*(mhist.^3).*hhist; 
  gk=gkmax*nhist.^4;
  clf
  plot(t,gna,'r');
  hold on
  plot(t,gk,'b');
  legend('gna','gk')
  hold off
  
figure(4);
X=0:i_step_size:0.6;
plot(X,no_of_peaks*1000/(niter/100));
xlabel('I Applied');
ylabel('Freq of Spikes')
hold on
for l=2:length(no_of_peaks) 
    if no_of_peaks(l)>0 && no_of_peaks(l-1)==0       
        I1=(l-1)*i_step_size;
    end
    if no_of_peaks(l)> no_of_peaks(l-1)+4        
        I2=(l-1)*i_step_size;
    end
    if no_of_peaks(l)< no_of_peaks(l-1)-  2        
        I3=(l-2)*i_step_size;
    end
end
I1=0*no_of_peaks*1000/(niter/100)+I1;
plot(I1,no_of_peaks*1000/(niter/100),'r');  
text(I1(1),-3,'I1');
I2=0*no_of_peaks*1000/(niter/100)+I2;
plot(I2,no_of_peaks*1000/(niter/100),'g');
text(I2(1),-3,'I2');
I3=0*no_of_peaks*1000/(niter/100)+I3;
plot(I3,no_of_peaks*1000/(niter/100),'y'); 
text(I3(1),-3,'I3');
text(0.5,100,'Before I1 => No A.P.s are seen');
text(0.5,90,'I1 to I2 => Finite no. of  A.P.s are seen');
text(0.5,80,'I2 to I3 => Infinite no.of A.P.s are seen');
text(0.5,70,'After I3 => No A.P.s are seen');
disp(I1(1))
disp(I2(1))
disp(I3(1))