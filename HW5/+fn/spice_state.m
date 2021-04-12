function rv = spice_state(epoch, target, frame, abcorr, observer) 

    rv = zeros(length(epoch), 6); 
    
    for i = 1:length(epoch) 

        %  Look-up the state for the defined parameters.
        starg   = mice_spkezr( target, epoch(i), frame, abcorr, observer);
        rv(i,:) = starg.state(1:6); 
        
    end 

end 