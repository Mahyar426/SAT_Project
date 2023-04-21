% Sample code for NMS iterative decoding

clc
clear
close all
load 128_64_LDPCcode.mat

for row = 1:size(H,1)
    Tanner_v2c{row} = find(H(row, :));
end

for col = 1:size(H,2)
    Tanner_c2v{col} = (find(H(:, col)))';
end

Out = (In < 0);
sindrome = mod(H * Out, 2);   

if sum(sindrome)~=0    
          
    [m, n] = size(H);
    Post=zeros(n,1);
    NMSfactor = 0.8;
           
    LLR = In;   % LLR per MS e NMS
    %LLR = 2*In/var;   % LLR per SPA
        
    Mess=H.*LLR';    
                        
    for Iter=1:IterMax            
        for check = 1:m         
            variab_to_check_mess = Mess(check, Tanner_v2c{check});
            for t=1:length(variab_to_check_mess)
                segni        = sign(variab_to_check_mess);
                magnitude        = abs(variab_to_check_mess);
                segni(t)     = 1;
                magnitude(t)     = inf;
                check_to_variab_mess(t) = prod(segni) * min(magnitude) * NMSfactor;
            end %t

            Mess(check, Tanner_v2c{check}) = check_to_variab_mess;    
    
        end % check

       

        for variab=1:n  
                                 
            check_to_variab_mess = Mess(Tanner_c2v{variab},variab);
            
            Poste = LLR(variab) + sum(check_to_variab_mess);
            variab_to_check_mess = (Poste - check_to_variab_mess)*NMSfactor;           
            Mess(Tanner_c2v{variab}, variab)   = variab_to_check_mess;   
            Post(variab)=Poste;     
        end % variab

                
        Out     = double(Post < 0);

        sindrome              = mod(H * Out, 2);
    
        if sum(sindrome)==0
              break;
        end
    end % for num_iter   
end % if syndrome not 0