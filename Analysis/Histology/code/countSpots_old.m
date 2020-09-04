function countSpots_old(imPath, wellPath, sgeID)

  d = dir(fullfile(imPath, '*.mat'));
  fname = d(sgeID).name;
  load(fullfile(imPath, fname));
  x = strsplit(fname, '_');
  jsonname = fullfile(wellPath, x{1}, 'scalefactors_json.json');
  w = jsondecode(fileread(jsonname));
  R = ceil(w.spot_diameter_fullres/2);
  tbl = readtable(fullfile(wellPath, x{1}, 'tissue_positions_list.txt'));
  nSpots = size(tbl, 1);
  count = zeros(nSpots, 1);
  mask = zeros(size(BW));
  crow = table2array(tbl(:, 5));
  ccol = table2array(tbl(:, 6));
  for i = 1:nSpots
    mask(crow(i), ccol(i)) = 1;
  end
  mask = bwdist(mask) <= R;
  mask = bwlabel(mask);
  parfor i = 1:nSpots
      disp(i)
      idx = mask(crow(i), ccol(i));
      tmpBW = BW;
      tmpBW(mask~=idx) = 0;
      [~, c] = bwlabel(tmpBW);
      count(i) = c;
  end
  tbl = [tbl array2table(count)];
  tbl.Properties.VariableNames = {'barcode','tissue','row','col','imagerow','imagecol','count'};
  mkdir(fullfile('/users/jcatalli/Histology', x{1}));
  writetable(tbl, fullfile('/users/jcatalli/Histology', x{1}, 'tissue_spot_counts1.csv'), 'Delimiter', ',');

end
