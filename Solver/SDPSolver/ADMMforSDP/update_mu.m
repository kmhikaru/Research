function [myu,it_pinf,it_dinf] = update_mu(pinf,dinf,eta1,eta2,it_pinf,it_dinf,h4,gmm,myu,myu_min,myu_max)
  if pinf/dinf <= eta1
    it_pinf = it_pinf + 1;
    it_dinf = 0;
    if it_pinf >= h4
      myu = max([gmm * myu,myu_min]);
      it_pinf = 0;
    end
  end
  if pinf/dinf > eta2
    it_dinf = it_dinf + 1;
    it_pinf = 0;
    if it_dinf >= h4
      myu = min([myu/gmm,myu_max]);
      it_dinf = 0;
    end
  end
