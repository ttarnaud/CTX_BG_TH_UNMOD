function CONN = SynC(N_pre,N_post,C)
% C is a non-zero integer, indicating how much synapses converge on a
% single post-synaptic neuron
C = min(C,N_pre);
CONN = gallery('circul',[ones(C,1); zeros(N_pre-C,1)])';
if N_pre>=N_post
CONN = CONN(:,1:N_post);
else
CONN = repmat(CONN,1,ceil(N_post/N_pre));
CONN = CONN(:,1:N_post);
end
end