function [t_sample, loopedAround] = incrementTPattern(t_sample, nClasses)
  pos = 1;
  loopedAround = false;
  
  while(pos>0)
    t_sample(pos) = t_sample(pos)+1;
    if t_sample(pos) > nClasses
      t_sample(pos) = 1;
      if numel(t_sample)>pos
        pos = pos+1;
      else
          t_sample = 1;
          loopedAround = true;
      end
    else
      pos = 0;
    end
  end
end