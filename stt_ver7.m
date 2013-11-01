function res=stt_ver7(n,i)

% Choice 1
seq=rand(n,1);
seq(find(seq<0.3))=0;
seq(find(seq>=0.3))=1;
%seq=[1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1];
seq1=zeros(1,150);
seq2=zeros(1,25);
seq=[seq1 seq2 seq1 seq2 seq1 seq2 seq1 seq2];
aux=sum(seq==i)/length(seq);
res1=(1/(1-aux))

% Choice 2
aux2=0;
aux3=0;
count=0;
for j=2:length(seq)
    if ((seq(j)==i) && (seq(j)==seq(j-1)))
        aux2=aux2+1;
    else if (((seq(j-1)==i) && ~(seq(j)==i)))
            aux3=aux3+aux2+1;
            count=count+1;
            aux2=0;
        end
    end
    if (seq(j)==i) && (j==length(seq))
        aux3=aux3+aux2+1;
        count=count+1;
    end
end

res2=aux3/count

end