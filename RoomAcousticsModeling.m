function [ newCoord, nextside ] = nextIntersect(previousIntersect, vector, plane, normal )
% this function finds the ray next intersection given the ray previous 
%intersection, current side, ray vector,and planes matrix
intersect = zeros(3,6);  state = zeros(6,3);  Dn = zeros(6,1); newCoord = zeros(3,1);t = zeros(6,1);
normal2 = -1*normal ; % to get the outside normal vector 
nextside=7;
vector =roundn(vector, -6);
for i =1:6
    Dn(i) = normal2(i,:)*vector;
    if Dn(i) > 0
        t(i,1) = ((plane(i,:) - previousIntersect')*normal2(i,:)')/Dn(i); 
        intersect(:,i) = t(i)*vector+ previousIntersect;
        intersect(:,i) = roundn(intersect(:,i),-6);
        for k = 1:3
            if intersect(k,i) >= min(plane(:,k)) && intersect(k,i) <= max(plane(:,k)) 
                state(i,k) = 1;                
            end
        end
    end
end

idx1 = Dn > 0;
idx2 = t >= 0 ;
idx = idx1.*idx2;

for i =1 :6 
    if sum(state(i,:)) == 3 && t(i) == min(t(logical(idx)))
        nextside = i; 
        newCoord = intersect(:,i);
        return
    end   
end
end

x = 224; % relative number of phi and tetha angles
t=linspace(0,180,x); p = linspace(0,359,2*x);
maxReflection =70; % the order of reflection
[T, P ]= meshgrid(t, p);% for each value of tetha, phi will change all the way fr om 0 to 359
T=T(:)'; P= P(:)'; 
N = length(t)*length(p); % total number of rays
r = zeros(3, N, maxReflection );
r(:,:,1)=[sind(T).*cosd(P); sind(T).*sind(P); cosd(T)];% create unit vector rays
speed = 343; dist = zeros(maxReflection,N); % distance traveld in each step and the total distance
d = zeros(maxReflection+1,N); % to indicate which side the ray will land on 1-6
coord = zeros(size(r)); % the (x,y,z) oordinates of the orgin and the distination of each ray 

normal = [-1 0 0 ; 0 -1 0; 1 0 0 ; 0 1 0 ; 0 0 1; 0 0 -1]; % the rows are normal vector for each side [1 - 6] normal inside

%% room dimension
L = 12; W = 10; H = 3; numReceivers = 3;
S = [2;5;2]; R = [6 5 1.5; 10 5 1.5; 11 9 1.5]'; R=R-S*ones(1,numReceivers);  
f = [125 	250 	500 	1000 	2000 	4000 ]';
a2 =  -S(1); a1 = L -S(1); 
b2 = -S(2); b1 = W -S(2); 
c2 = -S(3); c1 = H -S(3); 
vRoom = L*W*H;
vReceiver = 10*vRoom/N;
radius = (vReceiver*0.75/pi)^(1/3) ;

lMax = sqrt(L^2 + W^2 + H^2);
received = zeros(maxReflection,N,numReceivers,length(f));
timeMatrix = zeros(maxReflection,N,numReceivers);
absorptionCoff = [0.22 0.2 0.2 0.2 0.15 0.15; 0.22 0.2 0.2 0.2 0.15 0.15; 0.22 0.2 0.2 0.2 0.15 0.15; 0.22 0.2 0.2 0.2 0.15 0.15; 
    0.15 0.15 0.1 0.1 0.1 0.05; 0.2 0.2 0.22 0.2 0.2 0.1];  

refCoff = (1 -absorptionCoff ).^0.5; %reflection coffecients for each side [1-6] for differnt freqs
plane = [a1 0 0; 0 b1 0; a2 0 0; 0 b2 0 ; 0 0 c2; 0 0 c1]; % each row defined the planes corresponding to the sides of the room [1 - 6]

air_coff = [0.445	1.32	2.73	4.66	9.86	29.4]; %125 Hz	250 Hz	500 Hz	1000 Hz	2000 Hz	4000 Hz	db per Km
air_coff = air_coff/1000;
rayEnergy = zeros( maxReflection,N,numReceivers,length(f));% rayEnergy(1,:) = E0;
rayLevel = zeros(maxReflection,N,length(f)); rayLevel(1,:,:) = 1; % array for the rays signal levels
totalDist = zeros(N,1);
Lw = 90;   E = zeros(maxReflection, N); W0 = (10^(Lw/10)/N)*1e-12;
scale = W0/vReceiver;
dist_inter = zeros(maxReflection, N);

energyLimit = zeros(length(f),1);
energyValue = zeros(N,length(f));
Si = [W*H L*H W*H L*H L*W L*W];
vRoom = L*W*H;
den = zeros(length(f),1);
den2 = zeros(length(f),1);
T60 = zeros(length(f),1);
t60 = zeros(length(f),1);

for i = 1: length(f)
        den(i) = Si*absorptionCoff(:,i);
        T60(i) = 0.161*vRoom/den(i);
        den2(i) = -Si*log(1 - absorptionCoff(:,i));
        t60(i) = 0.161*vRoom/(den2(i) + 4*air_coff(i)*vRoom);
        energyLimit(i) = (scale*abs(lMax*air_coff(i)))/1e6;
end

scaleDot = zeros(numReceivers,1);
for k = 1: numReceivers
    scaleDot(k,1) = R(:,k)'*R(:,k);
end

%% for room visualization 
res = 100; scaleZ = H; scaleY =2;
range_x = linspace(a2,a1,res)'; X1 = a1*ones(res,1); X2 = a2*ones(res,1); Y2 = b2*ones(res,1); Y1 = b1*ones(res,1); Z1= c1*ones(res,1); Z2= c2*ones(res,1);
walls = [a2 a1; b2 b1; c2 c1]; range = zeros(3,res);
for i = 1: 3
   range(i,:) = linspace(walls(i,1),walls(i,2),res);
end
wallNormals = [-1 0 0; 0 -1 0 ; 1 0 0; 0 1 0; 0 0 -1; 0 0 1];
quiver3(X1, range(2,:)', Z2,zeros(res,1), zeros(res,1), scaleZ*Z1,':', 'AutoScale','off','ShowArrowHead','off');
hold on
quiver3(range(1,:)', Y1, Z2,zeros(res,1), zeros(res,1),scaleZ* Z1,':', 'AutoScale','off','ShowArrowHead','off');
quiver3(X2, range(2,:)', Z2,zeros(res,1), zeros(res,1),scaleZ* Z1,':', 'AutoScale','off','ShowArrowHead','off');
quiver3(range(1,:)', Y2, Z2,zeros(res,1), zeros(res,1),scaleZ* Z1,':', 'AutoScale','off','ShowArrowHead','off');
quiver3(range(1,:)', Y1, Z1,zeros(res,1), scaleY*Y2,zeros(res,1),':', 'AutoScale','off','ShowArrowHead','off');
quiver3(range(1,:)', Y1, Z2,zeros(res,1), scaleY*Y2,zeros(res,1),':', 'AutoScale','off','ShowArrowHead','off');
scatter3(0,0,0,140,'y','filled')
scatter3(R(1,:),R(2,:),R(3,:),[140 140 140 ]','m','filled')
text(0, 0,0, 'S')
for i =1: 3
    text(R(1,i), R(2,i), R(3,i), sprintf('R%d', i))
end
xlabel('Length = 12 [m]')
ylabel('Width = 10 [m] ')
zlabel('Height = 3 [m]')
title('Room 3D view')


%%
% for the intial emission of rays
for i =1: N  
    [ coord(:, i, 2), d(1,i)] = nextIntersect([ 0 0 0 ]', r(:,i,1), plane, normal); 
    dist(1,i) = norm(coord(:,i,2));

    for k =1:numReceivers
        temp1 = r(:,i,1)'*(coord(:,i,1) - R(:,k));
        temp = temp1^2 - (scaleDot(k,1) + coord(:,i,1)'*coord(:,i,1) -2*coord(:,i,1)'*R(:,k) - radius^2);

        if temp > 0 && (R(:,k)-coord(:,i,1))'*r(:,i,1) > -(R(:,k)-coord(:,i,1))'*r(:,i,1) 
             scale1 = -temp1 - sqrt(temp);
             scale2 = -temp1 + sqrt(temp);
             coord0 = coord(:,i,1) + scale1*r(:,i,1); 
             coord1 = coord(:,i,1) + scale2*r(:,i,1); 
             intersectLength = norm(coord1 - coord0);
             for l = 1: length(f)
                received(1, i,k,l ) = 1;
                timeMatrix(1,i,k) = (norm(R(:,k)))/speed;
                rayEnergy(1,i,k,l) = scale*rayLevel(1,i,l)*intersectLength*exp(-norm(coord(:,i,1) - R(:,k))*air_coff(l));
             end
        end
    end     
end

for l =1: length(f)
    energyValue(:,l) = radius*scale*rayLevel(1,:,l).*exp(-dist(1,:)*air_coff(l));  
end

%% for the subsequent reflection of rays
state = zeros(3,1);
for j = 2:maxReflection
  for l = 1: length(f)
    for i =1 :N 
        if abs(energyValue(i,l)) > energyLimit(l,1)
        rayLevel(j,i,l) = rayLevel(j-1,i,l);
        state(:) = 0;
         
           if abs(coord(1,i,j)- a1) < 1e-6 || abs(coord(1,i,j)- a2) < 1e-6
               r(:,i,j) = r(:,i,j-1).*[-1;1;1];
               rayLevel(j,i,l) = rayLevel(j,i,l)*((coord(1,i,j)==a1)*(1-absorptionCoff(1,l)) + (coord(1,i,j)==a2)*(1-absorptionCoff(3,l)));
               state(1,1) = 1;
           end

           if abs(coord(2,i,j)- b1) <  1e-6 || abs(coord(2,i,j)- b2) < 1e-6 
               r(:,i,j) = r(:,i,j-1).*[(-1)^state(1,1);-1;1];
               rayLevel(j,i,l) = rayLevel(j,i,l)*((coord(2,i,j)==b1)*(1-absorptionCoff(2,l)) + (coord(2,i,j)==b2)*(1-absorptionCoff(4,l))); 
               state(2,1) = 1;
           end

           if abs(coord(3,i,j)- c1) < 1e-6  || abs(coord(3,i,j)- c2) < 1e-6 
               r(:,i,j) = r(:,i,j-1).*[(-1)^state(1,1);(-1)^state(2,1);-1];
               rayLevel(j,i,l) = rayLevel(j,i,l)*((coord(3,i,j)==c1)*(1-absorptionCoff(6,l)) + (coord(3,i,j)==c2)*(1-absorptionCoff(5,l)));
               state(3,1) = 1;
           end       
           [ coord(:, i, j+1), ~] = nextIntersect(coord(:, i, j ), r(:,i,j), plane, normal);  
           rayLevel(j,i,l) = abs(rayLevel(j,i,l))*(-1)^(j-1);
           dist(j,i) = norm(coord(:,i,j+1) - coord(:,i,j)) + dist(j-1,i);           
           for k =1:numReceivers
                temp1 = r(:,i,j)'*(coord(:,i,j) - R(:,k));
                temp = temp1^2 - (scaleDot(k,1) + coord(:,i,j)'*coord(:,i,j) -2*coord(:,i,j)'*R(:,k) - radius^2);
                if temp > 0 && (R(:,k)-coord(:,i,j))'*r(:,i,j) > -(R(:,k)-coord(:,i,j))'*r(:,i,j) 
                    scale1 = -temp1 - sqrt(temp);
                    scale2 = -temp1 + sqrt(temp);
                    received(j, i,k,l ) = 1;
                    coord0 = coord(:,i,j) + scale1*r(:,i,j); 
                    coord1 = coord(:,i,j) + scale2*r(:,i,j); 
                    intersectLength = norm(coord1 - coord0);
                    timeMatrix(j,i,k,l) = (dist(j-1,i) + norm(coord(:,i,j) - R(:,k)))/speed;
                    rayEnergy(j,i,k,l) = scale*abs(rayLevel(j,i,l))*intersectLength*exp(-(dist(j-1,i) + norm(coord(:,i,j) - R(:,k)))*air_coff(l));
                end
           end
       end  
    end
    energyValue(:,l) = radius*scale*rayLevel(j,:,l).*exp(-dist(j,:)*air_coff(l));  
  end

end

%% for saving and loading the parameters
%save('variabls')
%load('variabls')

%%
fs = 44100; % sampling frequency
maxTimeConsidered = 1;
deltaTime = 1/fs; timeIndex = 0:deltaTime: maxTimeConsidered +deltaTime;
timeIndex = timeIndex';

newEnergy = zeros(numReceivers,length(timeIndex),length(f)); amplitude =zeros(numReceivers,length(timeIndex),length(f));
newEnergy2 = zeros(numReceivers,length(timeIndex),length(f));
signalSign =zeros(numReceivers,length(timeIndex),length(f));
% for obtaining the RIR
for i = 1: numReceivers
    for l = 1:length(f)
    id1=(received(:,:,i,l)==1); 
    tt = timeMatrix(:,:,i,l);  
    receivedSignal.(sprintf('time%d',i)) = tt(logical(id1));
    ray1 = rayLevel(:,:,l);
    ray = ray1(logical(id1));
    energy1 = rayEnergy(:,:,i,l);
    energy =energy1(logical(id1));
    for j = 1: length(timeIndex)-1
        tempIndex1 = receivedSignal.(sprintf('time%d',i)) > timeIndex(j);
        tempIndex2 = receivedSignal.(sprintf('time%d',i)) <= timeIndex(j+1);
        tempIndex = tempIndex1.*tempIndex2;
        tempVector = (energy(logical(tempIndex)));
        tempVector2 = (ray(logical(tempIndex)));
    
         if isscalar(max(tempVector)) == 1  
             [newEnergy(i,j,l)] = sum(abs(tempVector));
             signalSign(i,j,l) = sign(max(tempVector2));      
         end
    end
    tempMax = max(newEnergy(i,:,l));
    energyFlag1 = abs(newEnergy(i,:,l)) < tempMax/1e6 ;
    energyFlag2 = abs(newEnergy(i,:,l)) > 0;
    energyFlag = energyFlag1.*energyFlag2;
    [~, endIndex] = find(energyFlag>0);
    if isscalar(endIndex) == 0
        endIndex = length(newEnergy);
    end
    newEnergy2(i,1:endIndex(1),l) = newEnergy(i,1:endIndex(1),l);
    amplitude(i,:,l) = signalSign(i,:,l).*sqrt(400*abs(newEnergy2(i,:,l)));
    amplitude(i,:,l) = amplitude(i,:,l)/max(abs(amplitude(i,:,l)));
    end
end
%% for estimating the SPL 
SPLarray = zeros(numReceivers, length(f));
for i = 1: numReceivers
    for l = 1:length(f)
        SPLarray(i,l) = 10*(log10(sum(newEnergy2(i,:,l)))+12);        
    end
    hold on
    scatter(f,SPLarray(i,:),'filled')
end
legend('R1','R2','R3')
title('SPL for each receiver at different frequencies')
ylabel('SPL [dB]')
xlabel('Frequency [Hz]')
grid on
grid minor
ylim([65 80])
xlim([0 4500])
receiverSPL = sum(SPLarray')/length(f);
hold off

%% manual calculation for the SPL at each receiver at each frequency band
S = 2*(L*W + H*(W+L)) ; Q =4; alpha = zeros(numReceivers,1);
Lp = zeros(numReceivers,1);
 
for i =1: numReceivers
    alpha(i) = sum(absorptionCoff(i,:))/length(f);
    Lp(i) = Lw + 10*log10(Q/(4*pi*norm(R(:,i))^2) + 4*(1-alpha(i))/(S*alpha(i)));  
end


scatter((1:numReceivers)',receiverSPL','filled')
hold on
scatter((1:numReceivers)',Lp','filled')
for i =1 : numReceivers
    text(i,receiverSPL(i)-1, ['delta=' num2str(Lp(i)-receiverSPL(i) ) 'dB']);
end
legend('Modelled SPL','Calculated SPL')
title('SPL at each receiver')
ylabel('SPL [dB]')
xlabel('Receiver ')
grid on
grid minor
ylim([65 80])
xlim([0 4])

%% for visualizing the impulse response for all receivers
freq = 1;
stopIndex = 25000;
subplot(3,1,1)
plot(timeIndex(1:stopIndex,1),(amplitude(1,1:stopIndex,freq)))
title('Impulse response from each receiver for 1KHz frequency')
ylabel('Relative amplitude R1')
xlabel('Time [seconds]')
subplot(3,1,2)
plot(timeIndex(1:stopIndex,1),(amplitude(2,1:stopIndex,freq)))
ylabel('Relative amplitude R2')
xlabel('Time [seconds]')
subplot(3,1,3)

plot(timeIndex(1:stopIndex,1),(amplitude(3,1:stopIndex,freq)))
ylabel('Relative amplitude R3')
xlabel('Time [seconds]')

%% for estimating the RT60
RT60 = zeros(numReceivers,length(f));

figure(1)
for i = 1: 3 %numReceivers
    for j = 1:length(f)
        figure(j)
        ind2 = newEnergy(i,:,j) ~= 0; 
        decayTime = timeIndex(ind2,1)';
        Idb = 10*log10(abs(newEnergy(i,ind2,j))/1e-12);
        decayCurve = Idb(1,:)-max(Idb(1,:));
        n=length(decayCurve);
        H=fftshift(fft(decayCurve,n));
        H([1: floor(n*0.49), floor(n*0.51):end]) = 0;
        h_rec=(ifft(H,n));
        [maxValue,startIndex] = max(-abs(h_rec));
        newCurve = -abs(h_rec(1,:))-maxValue;
        shiftUp=10; 
        smo = fit(decayTime(1,startIndex:end)',newCurve(1,startIndex:end)'+shiftUp, 'poly1');
        plot(smo,decayTime(1,startIndex:end-100)',newCurve(1,startIndex:end-100)');
        ylabel('Relative SIL [dB]')
        xlabel('Time [seconds]')
        title(sprintf('RT60 Receiver %d  Frequency %d Hz',i,f(j)))
        saveas(gcf,sprintf('RT60_R %d_F %d.jpg',i,j));
        fittedCurve = smo(decayTime(1,startIndex:end)');
        [time, ~] = find(fittedCurve <= -30);
        RT60(i,j) = 2*(decayTime(1,time(1)) -decayTime(1,startIndex));
        
    end
end
%%
meanT60 = mean(RT60);
100*(abs(meanT60 -t60')./t60')
errorbar(f,t60',t60'-meanT60)
scatter(f, t60')
hold on
scatter(f, meanT60,'*')
legend('Calculated','Mean modelled T60')
ylim([0.5  1])
xlim([0 4500])
ylabel('T60 [seconds]')
xlabel('Frequency [Hz]')

%% for combining the frequency response from different octave bands

fftSize = length(timeIndex);length(timeIndex);freq=1; Receiver=1;
octaveBW = [f(1:end-1)*sqrt(2)];
combinedFR = zeros(numReceivers, fftSize/2);
freqRange = 0:fs/fftSize:fs/2 -1;
FR_BW = zeros(length(octaveBW),1);
for j = 1:length(octaveBW)
    temp = freqRange >= octaveBW(j);
    [~, ind] = find(temp ==1);
    FR_BW(j,1) = ind(1);
end

FR_BW = [1;FR_BW;length(freqRange)];
orignalFFT = zeros(numReceivers, fftSize);

for i = 1: numReceivers
    for j = 1:length(f)
        Hdft = fftshift(fft((amplitude(i,:,l)),fftSize));        
        HdftdB = (20*log10(abs(Hdft)));
        HdftdBF = fftshift(fft(HdftdB,length(HdftdB)));
        HdftdBsmooth =HdftdBF;
        HdftdBsmooth([1: floor(fftSize*0.495), floor(fftSize*0.505):end])=0;
        smoothFR = ifft(HdftdBsmooth,fftSize);
        
        for k = 1:length(FR_BW)-1        
            combinedFR(i,FR_BW(k):FR_BW(k+1))= smoothFR([floor(fftSize*0.5)+FR_BW(k):floor(fftSize*0.5)+FR_BW(k+1)]);
            orignalFFT(i,FR_BW(k):FR_BW(k+1)) = Hdft([floor(fftSize*0.5)+FR_BW(k):floor(fftSize*0.5)+FR_BW(k+1)]);
        end
    end
end

combindFFT = [fliplr(combinedFR(:,2:end)) combinedFR];
orignalFFT = [fliplr(orignalFFT(:,2:end)) orignalFFT];
semilogx(freqRange,abs(combinedFR(1,1:end-1))-max(abs(combinedFR(1,1:end-1))))

subplot(3,1,1)
semilogx(freqRange,abs(combinedFR(1,1:end-1))-max(abs(combinedFR(1,1:end-1))))
title('Relative frequency response for each receiver ')
ylabel('Freq. Response R1 [dB]')
xlabel('Frequency [kHz]')
ylim([-30,5])
xlim([0, 25000])

subplot(3,1,2)
semilogx(freqRange,abs(combinedFR(2,1:end-1))-max(abs(combinedFR(2,1:end-1))))
ylabel('Freq. Response R2 [dB]')
xlabel('Frequency [kHz]')
ylim([-30,5])
xlim([0, 25000])

subplot(3,1,3)
semilogx(freqRange,abs(combinedFR(3,1:end-1))-max(abs(combinedFR(3,1:end-1))))
ylabel('Freq. Response R3 [dB]')
xlabel('Frequency [kHz]')
ylim([-30,5])
xlim([0, 25000])

%% to get the imnpulse response corresponding to the combined FR for each receiver
finalIR = zeros(numReceivers,fftSize);
for i =1:numReceivers
    finalIR(i,:) = (ifft(orignalFFT(i,:),fftSize));  
    finalIR(i,:) =  finalIR(i,:)/max(abs(finalIR(i,:)));
end
finalIR = -real(fliplr(finalIR));
plot(timeIndex(1:fftSize,1),(finalIR(3,:)))

subplot(3,1,1)
plot(timeIndex(1:fftSize,1),(finalIR(1,:)))
title('Impulse response from each receiver')
ylabel('Relative amplitude R1')
xlabel('Time [seconds]')
subplot(3,1,2)
plot(timeIndex(1:fftSize,1),(finalIR(2,:)))
ylabel('Relative amplitude R2')
xlabel('Time [seconds]')
subplot(3,1,3)
plot(timeIndex(1:fftSize,1),(finalIR(3,:)))
ylabel('Relative amplitude R3')
xlabel('Time [seconds]')

%%
finalIR2= zeros(size(finalIR));
for i =1: numReceivers
        tempIndex = abs(finalIR(i,:)) > 1e-3;
        finalIR2(i,tempIndex) = finalIR(i,tempIndex);
end
%% for sound reproduction inside the room
receiver = 1;
figure
res = abs(finalIR(receiver,:))/max(abs(finalIR(receiver,:)));
SPLdB = 20*log10(res);
plot(timeIndex,SPLdB)
[s, fs]=audioread('dryspeech.wav');
h=finalIR2(receiver,:)';
h=h/max(h);
h=h(1:17000,1);
N=length(s)+length(h)-1;
window = hamming(2*length(h));
window2= window((length(h)):end-1);
h = h.*window2 ;
s= s.*(hamming(length(s)));
sound(s,fs)
pause(5)
figure(1)
h = [h; zeros(length(s)-1,1)]; s = [s ; zeros(length(h)-1,1) ]; 
srDFT = (fft(s(:,1), N));
hr = (fft(h(:,1),N));
output_r = srDFT.*hr; 
res_r = ifft(output_r); 
res = [res_r res_r ]; 
res = res/max(abs(res(:)));
freq = linspace(-fs/2, fs/2, N)';
sound(res,fs)
audiowrite(sprintf('ReproductionR%d.wav',receiver), res, fs)








