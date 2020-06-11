%     ###### #     #  ####      ####  #  ####  ####  ###  #     #                   ###    
%     #      #     # #    #     #   # # #     #     #   # #     #                  #   #   
%     #      #     # #    #     #   # # #     #     #   # #     #                  #   #   
%     #####  #     # ######     ####  # #     #     #   # #     #                   ###   #
%     #      #     # #    #     #     # #     #     #   # #     #                  #   # # 
%     #      #     # #    #     #     # #     #     #   # #     #                  #    ## 
%     ###### ##### # #    #     #     #  ####  ####  ###  ##### #                   ####  #


%     #     #  ####  ####  #  ####  #   #      #### #####  ####  #####  ####   #### #   # #####
%     ##   ## #    # #   # # #    # ##  #     #       #   #    #   #   #    # #     #   # #
%     # # # # #    # #   # # #    # # # #     #       #   #    #   #   #    # #     #   # #
%     #  #  # ###### ####  # ###### #  ##      ###    #   ######   #   ###### #     ##### ####
%     #     # #    # # #   # #    # #   #         #   #   #    #   #   #    # #     #   # #
%     #     # #    # #  #  # #    # #   #         #   #   #    #   #   #    # #     #   # #
%     #     # #    # #   # # #    # #   #     ####    #   #    #   #   #    #  #### #   # #####


% Anno Accademico 2019/2020
% Marian Statache, Elia Piccoli

% Script per il defect-detection nei tessuti. 

clear all;
close all;

list=dir('images/tex*');
s='images/';

for index=1:size(list,1)
    name=strcat(s,list(index).name);
    A=rgb2gray(imread(name));
    [R,C]=size(A);
    
    A_mod=abs(fft2(A));
    A_mod=(A_mod-min(A_mod(:)))/(max(A_mod(:))-min(A_mod(:)))*255;
    
    B=A_mod(1:R/2,1:C/2);
    B(1:10,:)=0;
    B(:,1:10)=0;
    [val, ind] = max(B(:));
    [RR CC]=size(B);
    [row col] = ind2sub([RR CC],ind);
    
    Xfactor=col/2;  %prova la metà del più piccolo di fft2
    Yfactor=row/2;
    
    if(Xfactor>32)
        Xfactor=32;
    elseif(Xfactor<4)
        Xfactor=4;
    end
    
    if(Yfactor>32)
        Yfactor=32;
    elseif(Yfactor<4)
        Yfactor=4;
    end
    
    % scelta dei pattern dinamica, sia per la dimensione che per la posizione
    
    pattern1 = A(1:R/Yfactor,1:C/Xfactor);
    pattern2 = A(3:R/Yfactor+2,3:C/Xfactor+2);
    pattern3 = A(R/(Yfactor*2):R/Yfactor+R/(Yfactor*2)-1,C/(Xfactor*2):C/Xfactor+C/(Xfactor*2)-1);
    pattern4 = A(1:R/Yfactor,C-C/Xfactor+1:C);
    pattern5 = A(3:R/Yfactor+2,C-C/Xfactor-1:C-2);
    pattern6 = A(R/(Yfactor*2):R/Yfactor+R/(Yfactor*2)-1,C-C/Xfactor+1-C/(Xfactor*2):C-C/(Xfactor*2));
    pattern7 = A(R-R/Yfactor+1:R,1:C/Xfactor);
    pattern8 = A(R-R/Yfactor-1:R-2,3:C/Xfactor+2);
    pattern9 = A(R-R/Yfactor+1-R/(Yfactor*2):R-R/(Yfactor*2),C/(Xfactor*2):C/Xfactor+C/(Xfactor*2)-1);
    
    
    c1 = normxcorr2(pattern1,A);
    c2 = normxcorr2(pattern2,A);
    c3 = normxcorr2(pattern3,A);
    c4 = normxcorr2(pattern4,A);
    c5 = normxcorr2(pattern5,A);
    c6 = normxcorr2(pattern6,A);
    c7 = normxcorr2(pattern7,A);
    c8 = normxcorr2(pattern8,A);
    c9 = normxcorr2(pattern9,A);
    
    c = (c1+c2+c3+c4+c5+c6+c7+c8+c9)/9;
    c = c(R/Yfactor:end-R/Yfactor+1,C/Xfactor:end-C/Xfactor+1);
    
    c=abs(c);
    
    % scelta della soglia in modo automatico
    soglia=mean(c(:));
    mask = c<soglia;

    % scelta automatica del disco
    if(soglia>0.1)
        disco=2;
    else
        disco=3;
    end
    
    se = strel('disk',disco);
    mask2 = imopen(mask,se);
    
    A=A(R/Yfactor-R/(Yfactor*2):end-R/Yfactor+R/(Yfactor*2),C/Xfactor-C/(Xfactor*2):end-C/Xfactor+C/(Xfactor*2));
    
    [RA CA] = size(A);
    [Rm Cm] = size(mask2);
    if(RA > Rm)
        A = A(1+(RA-Rm):end, :);
    else
        mask2 = mask2(1:end-(Rm-RA), :);
    end
    if(CA > Cm)
        A = A(:, 1+(CA-Cm):end);
    else
        mask2 = mask2(:, 1:end-(Cm-CA));
    end
    
    A1 = A;
    A1(mask2)=255;
    %le prossime 9 righe servono per rendere l'evidenziatura rossa anche
    %quando il difetto è bianco o tanto chiaro, in modo che si veda.
    %se si preferisce com'era prima basta scommentare la riga 108
    A2 = A;
    for i=1:size(A2,1)
        for j=1:size(A2,2)
            if(A2(i,j)>100&&A1(i,j)==255)
                A2(i,j)=100;
            end
        end
    end
    Af=cat(3,A1,A2,A2);
%     Af=cat(3,A1,A,A);
    figure(1);
    set(gcf, 'Position', get(0, 'Screensize'));
    subplot(5, 4, index);
    A_color=imread(name);
    A_color=A_color(R/Yfactor-R/(Yfactor*2):end-R/Yfactor+R/(Yfactor*2),C/Xfactor-C/(Xfactor*2):end-C/Xfactor+C/(Xfactor*2),:);
    imshowpair(A_color,Af,'montage')
    figure(2);
    imshowpair(A_color,Af,'montage')
    pause(2);
end

close(figure(2));

%%  test fft
% clear all;
% close all;
% 
% A = rgb2gray(imread('images/tex01.jpg'));
% [R,C]=size(A);
% A_mod=abs(fft2(A));
% A_mod=(A_mod-min(A_mod(:)))/(max(A_mod(:))-min(A_mod(:)))*255;
% 
% B=A_mod(1:R/2,1:C/2);
% B(1:6,:)=0;
% B(:,1:6)=0;
% imshow(B); title('TdF');
% [val, ind] = max(B(:));
% [R C]=size(B);
% [row col] = ind2sub([R C],ind);