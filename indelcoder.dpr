program indelcoder;

{$APPTYPE CONSOLE}

{$IFDEF FPC}
 {$MODE DELPHI}
{$ENDIF FPC}

uses
  SysUtils,Classes;

type
  Tgap = record
          s,e : integer;
        end;
  Tgaplist = array of Tgap;


function GetInputSequences(fname : string; seqs, taxa : TStringList) : boolean;
var
   infile : TextFile;
   instr : string;
   n,nt : integer;
   done,interleave : boolean;
begin
     result := true;
     interleave := false;
     AssignFile(infile,fname); Reset(infile);
     repeat
           readln(infile,instr);
     until eof(infile) or (pos('BEGIN TAXA',UpperCase(instr))>0)
           or (pos('BEGIN DATA',UpperCase(instr))>0);
     if eof(infile) then result := false
     else begin
          if (pos('BEGIN TAXA',UpperCase(instr)) > 0) then begin
             repeat
                   readln(infile,instr);
             until (Pos('NTAX',UpperCase(instr)) > 0);
             instr := Copy(instr,Pos('NTAX',UpperCase(instr))+5,length(instr));
             nt := StrToInt(Copy(instr,1,Pos(';',instr)-1));
             repeat
                   readln(infile,instr);
             until (Pos('TAXLABELS',UpperCase(instr)) > 0);
             instr := trim(Copy(trim(instr),Pos(' ',trim(instr)),length(instr)));
             for n := 1 to nt - 1 do begin
                 taxa.Add(trim(Copy(instr,1,Pos(' ',instr))));
                 instr := trim(Copy(instr,Pos(' ',instr),length(instr)));
             end;
             taxa.Add(Copy(instr,1,Pos(';',instr)-1));
             repeat
                   readln(infile,instr);
             until (Pos('MATRIX',UpperCase(instr)) > 0);
             for n := 1 to nt do begin
                 readln(infile,instr);
                 seqs.Add(trim(instr));
             end;
          end else begin
             repeat
                   readln(infile,instr);
             until (Pos('NTAX',UpperCase(instr)) > 0);
             instr := Copy(instr,Pos('NTAX',UpperCase(instr))+5,length(instr));
             nt := StrToInt(Copy(instr,1,Pos(' ',instr)-1));
             repeat
                   readln(infile,instr);
                   if (Pos('INTERLEAVE',UpperCase(instr)) > 0) then
                      interleave := true;
             until (Pos('MATRIX',UpperCase(instr)) > 0);
             for n := 1 to nt do begin
                 readln(infile,instr);
                 taxa.Add(Copy(instr,1,Pos(' ',instr)-1));
                 seqs.Add(trim(Copy(instr,Pos(' ',instr),length(instr))));
             end;
             if interleave then begin
                done := false;
                repeat
                      readln(infile,instr);
                      if (trim(instr) <> '') then begin
                         if (trim(instr) = ';') then done := true
                         else begin
                              seqs[0] := seqs[0] +
                                      trim(Copy(instr,Pos(' ',instr),length(instr)));
                              for n := 1 to nt - 1 do begin
                                  readln(infile,instr);
                                  seqs[n] := seqs[n] +
                                      trim(Copy(instr,Pos(' ',instr),length(instr)));
                              end;
                         end;
                      end;
                until done;
             end;
          end;
     end;
     CloseFile(infile);
end;

procedure WriteOutputFile(fname : string; seqs, taxa : TStringList;
          GapList : TGapList; n : integer);
var
   outfile : TextFile;
   i : integer;
begin
     AssignFile(outfile,fname); Rewrite(outfile);
     writeln(outfile,'#NEXUS');
     writeln(outfile,'BEGIN DATA;');
     writeln(outfile,'  DIMENSIONS NTAX='+IntToStr(taxa.count)+' NCHAR='
       +IntToStr(length(seqs[0]))+';');
     writeln(outfile,'  FORMAT MISSING=? DATATYPE=DNA GAP=- EQUATE="0=A 1=C";');
	   writeln(outfile,'  OPTIONS GAPMODE=MISSING;');
     writeln(outfile,'MATRIX');
     for i := 0 to taxa.Count - 1 do
         writeln(outfile,taxa[i] + '     ' + seqs[i]);
     writeln(outfile,';');
     writeln(outfile,'END;');
     writeln(outfile);
     writeln(outfile,'BEGIN SETS;');
     writeln(outfile,'  CHARSET Original=1-'+IntToStr(n)+';');
     writeln(outfile,'  CHARSET InDelChar='+IntToStr(n+1) + '-' +
        IntToStr(length(seqs[0]))+';');
     writeln(outfile,'END;');

     writeln(outfile);
     writeln(outfile,'[ Indel Character        Sites');
     writeln(outfile,'  ---------------   ---------------');
     for i := 0 to length(GapList) - 1 do begin
     writeln(outfile,'       '+format('%5d',[n+1+i])+'        ' +
         format('%6d',[GapList[i].s]) + ' - '+
         IntToStr(GapList[i].e));
     end;
     writeln(outfile,']');

     CloseFile(outfile);
end;

procedure CalcGapCodes(oseqs,nseqs : TStringList; DoSimple: boolean; var UniqueGaps : TGapList);
var
   i,g,ns,c,n : integer;
   tstr : string;
   SeqGaps : array of TGapList;
   IsNew : boolean;
begin
     // copy original sequences into new sequences
     nseqs.Clear;
     ns := oseqs.Count;
     for n := 0 to ns - 1 do nseqs.add(oseqs[n]);
     // mask trailing and leading gaps with ?
     for n := 0 to ns - 1 do begin
         tstr := nseqs[n];
         c := 1;
         repeat
               if (tstr[c] = '-') then tstr[c] := '?';
               inc(c);
         until (tstr[c] <> '-') or (c > length(tstr));
         c := length(tstr);
         repeat
               if (tstr[c] = '-') then tstr[c] := '?';
               dec(c);
         until (tstr[c] <> '-') or (c < 1);
         nseqs[n] := tstr;
     end;
     // find all gaps
     SetLength(SeqGaps,ns);
     for n := 0 to ns - 1 do begin
         SeqGaps[n] := nil;
         tstr := nseqs[n];
         c := 1;
         repeat
               if (tstr[c] = '-') then begin
                  g := length(SeqGaps[n]);
                  SetLength(SeqGaps[n],g+1);
                  SeqGaps[n,g].s := c;
                  repeat
                        inc(c);
                  until (tstr[c] <> '-');
                  SeqGaps[n,g].e := c - 1;
               end else inc(c);
         until (c = length(tstr));
     end;
     // find unique gaps
     UniqueGaps := nil;
     for n := 0 to ns - 1 do
         for g := 0 to length(SeqGaps[n]) - 1 do begin
             i := length(UniqueGaps);
             if (i = 0) then begin
                SetLength(UniqueGaps,i+1);
                UniqueGaps[i].s := SeqGaps[n,g].s;
                UniqueGaps[i].e := SeqGaps[n,g].e;
             end else begin
                 IsNew := true;
                 for c := 0 to i - 1 do
                     if (SeqGaps[n,g].s = UniqueGaps[c].s) and
                        (SeqGaps[n,g].e = UniqueGaps[c].e)
                     then IsNew := false;
                 if IsNew then begin
                    SetLength(UniqueGaps,i+1);
                    UniqueGaps[i].s := SeqGaps[n,g].s;
                    UniqueGaps[i].e := SeqGaps[n,g].e;
                 end;
             end;
         end;
     if DoSimple then begin // Simple Indel Coding
        // code the gaps
        for n := 0 to ns - 1 do
            for i := 0 to length(UniqueGaps) - 1 do begin
                IsNew := false; // use to identify matched gaps
                for g := 0 to length(SeqGaps[n]) - 1 do
                    if (UniqueGaps[i].s = SeqGaps[n,g].s) and
                       (UniqueGaps[i].e = SeqGaps[n,g].e)
                    then begin
                         nseqs[n] := nseqs[n] + '1';
                         IsNew := true;
                    end;
                if not IsNew then begin
                // The gap is not in the sequence; figure out if it is a superset
                //  of not
                  IsNew := false; // use to identify superset gaps
                  for g := 0 to length(SeqGaps[n]) - 1 do
                    if (UniqueGaps[i].s >= SeqGaps[n,g].s) and
                       (UniqueGaps[i].e <= SeqGaps[n,g].e)
                    then IsNew := true;
                  if IsNew then nseqs[n] := nseqs[n] + '-'
                  else nseqs[n] := nseqs[n] + '0';
                end;
            end;
     end else begin // Complex Indel Coding
         // analyze identified gaps
         //ComplexGapAnalysis(UniqueGaps);
         // code the gaps


     end;
     SeqGaps := nil;
end;

function AddIndelCode(infilen,outfilen : string; DoSimple : boolean) : boolean;
var
   oseqs,nseqs,taxa : TStringList;
   GapList : TGapList;
begin
     result := true;
     oseqs := TStringList.Create;
     nseqs := TStringList.Create;
     taxa := TStringList.Create;
     GapList := nil;
     if GetInputSequences(infilen,oseqs,taxa) then begin
        CalcGapCodes(oseqs,nseqs,DoSimple,GapList);
        WriteOutputFile(outfilen,nseqs,taxa,GapList,length(oseqs[0]));
     end else result := false;
     GapList := nil;
     oseqs.Free;
     nseqs.Free;
     taxa.Free;
end;


procedure WriteUsage;
begin
     writeln('Usage: indelcoder <input_file> <output_file>')
end;

var
   IsGood : boolean;
   infilen,outfilen : string;
begin
  writeln;
  writeln('IndelCoder');
  writeln('  Michael S. Rosenberg');
  writeln('  Copyright (c) 2005-2010');
  writeln('  Version 1.0.1.1');
  writeln;
  if (ParamCount <> 2) then WriteUsage
  else begin
       infilen := ParamStr(1);
       outfilen := ParamStr(2);
       IsGood := true;
       if (infilen = '') then begin
          IsGood := false;
          writeln('Error: No input file was specified.');
          WriteUsage;
       end else if (outfilen = '') then begin
          IsGood := false;
          writeln('Error: No output file was specified.');
          WriteUsage;
       end else if not FileExists(infilen) then begin
           IsGood := false;
           writeln('Error: "'+infilen+'" does not exist.');
       end;
       if IsGood then begin
          if not AddIndelCode(infilen,outfilen,true) then
             writeln('Error: Failed to code "'+infilen+'"')
          else writeln('Conversion complete');
       end;
  end;

end.
