function changed = patch_final_reshape(filename);
  fid = fopen(filename);
  % read first line
  line = fgetl(fid);
  C = {};
  changed = false;
  while ischar(line)
    % previous line
    line0 = line;
    % read next line
    line = fgetl(fid);
    if strcmp(line,'end')
      new_line0 = regexprep(line0,'reshape\((.*),([0-9]+),1\);','reshape($1,[],1);');
      if ~strcmp(line0,new_line0)
        changed = true;
        line0 = new_line0;
      end
    end
    C = {C{:} line0};
  end
  fclose(fid);

  if changed
    fid = fopen(filename,'w');
    fprintf(fid,'%s\n',C{:});
    fclose(fid);
  end
end
