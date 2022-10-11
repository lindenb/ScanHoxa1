/*
The MIT License (MIT)

Copyright (c) 2022 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.beust.jcommander.Parameter;
import com.github.lindenb.jvarkit.io.IOUtils;
import com.github.lindenb.jvarkit.lang.StringUtils;
import com.github.lindenb.jvarkit.samtools.SAMRecordDefaultFilter;
import com.github.lindenb.jvarkit.util.bio.AminoAcids;
import com.github.lindenb.jvarkit.util.bio.AcidNucleics;
import com.github.lindenb.jvarkit.util.bio.GeneticCode;
import com.github.lindenb.jvarkit.util.bio.SequenceDictionaryUtils;
import com.github.lindenb.jvarkit.util.bio.fasta.ContigNameConverter;
import com.github.lindenb.jvarkit.util.jcommander.Launcher;
import com.github.lindenb.jvarkit.util.log.Logger;
import com.github.lindenb.jvarkit.util.Counter;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.SequenceUtil;

public class ScanHoxa1 extends Launcher {
	private static final Logger LOG = Logger.build(ScanHoxa1.class).make();

	@Parameter(names={"--references"},description="path to references",required = true)
	private Path references;
	@Parameter(names={"--mapq"},description="mapq")
	private int mapq=10;
	@Parameter(names={"--min"},description="min occurence")
	private int min_occurence =0;

	@Override
	public int doWork(List<String> args) {
		final Map<Path,SAMSequenceDictionary> path2dict = new HashMap<>();
		final GeneticCode geneticCode = GeneticCode.getStandard();
		try {
			try(BufferedReader br= Files.newBufferedReader(references)) {
				String line;
				while((line=br.readLine())!=null) {
					LOG.info("dict for "+line);
					final Path faidx = Paths.get(line);
					final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(faidx);
					path2dict.put(faidx, dict);
					}
				}
			
			String referenceSeq = null;
			System.out.println("#sample\tbamPath\tinsertion_count\tdeletion_count\tbase_at_27_135_314\thDNA\tpep\tsymbol\tseen.X.times\tfreq");
			
			for(final Path bamPath:IOUtils.unrollPaths(args)) {
				//if(bamPath.toString().contains("00PUM")==false) continue;
				LOG.info(bamPath);
				final SAMSequenceDictionary dict = SequenceDictionaryUtils.extractRequired(bamPath);
				Path faix=null;
				for(Map.Entry<Path, SAMSequenceDictionary> pair: path2dict.entrySet()) {
					if(SequenceUtil.areSequenceDictionariesEqual(dict,pair.getValue())) {
						faix = pair.getKey();
						break;
						}
					}
				if(faix==null) {
					LOG.info("cannot find ref for " + bamPath);
					continue;
					}
				final ContigNameConverter convert = ContigNameConverter.fromOneDictionary(dict);
				final String contig = convert.apply("chr7");
				final int chromStart = 27_135_310;
				final int chromEnd = 27_135_339;
				final int C_arg_pos = 27_135_314;
				if(StringUtils.isBlank(contig)) continue;
				if(referenceSeq==null) {
					try(ReferenceSequenceFile f=ReferenceSequenceFileFactory.getReferenceSequenceFile(faix)) {
						referenceSeq= f.getSubsequenceAt(contig, chromStart, chromEnd).getBaseString().toUpperCase();
					}
				}
				
				final SamReaderFactory srf = super.createSamReaderFactory();
				srf.referenceSequence(faix);
				try(SamReader sr = srf.open(bamPath)) {
					final Counter<String> counter = new Counter<>();
					final String sample = sr.getFileHeader().getReadGroups().
							stream().map(S->S.getSample()).
							filter(S->!StringUtils.isBlank(S)).
							findFirst().orElse(IOUtils.getFilenameWithoutCommonSuffixes(bamPath));
					try(CloseableIterator<SAMRecord> iter=sr.query(contig,chromStart,chromEnd,false)) {
						while(iter.hasNext()) {
							final SAMRecord rec = iter.next();
							if(!SAMRecordDefaultFilter.accept(rec,this.mapq)) continue;
							if(rec.getStart() > chromStart) continue;
							if(rec.getEnd() < chromEnd) continue;
							final Cigar cigar = rec.getCigar();
							final String bases = rec.getReadString();
							int refpos= rec.getAlignmentStart();
							int readpos=0;
							final StringBuilder sb=new StringBuilder();
							int insertion_count=0;
							int deletion_count=0;
							char base_at_27_135_314 = '.';
							for(CigarElement ce: cigar) {
								final CigarOperator op =ce.getOperator();
								final int len = ce.getLength();
								switch(op) {
									case P: break;
									case H: break;
									case S: readpos+=len;break;
									case N: case D:
										for(int i=0;i< len;i++) {
											if(refpos>=chromStart && refpos<=chromEnd) {
												deletion_count++;
												}
											refpos++;
											}
										break;
									case I:
										for(int i=0;i< len;i++) {
											if(refpos>=chromStart && refpos<=chromEnd) {
												sb.append(bases.charAt(readpos+i));
												insertion_count++;
												}
											}
										readpos+=len;
										break;
									case EQ: case M: case X:
										for(int i=0;i< len;i++) {
											if((refpos+i)>=chromStart && (refpos+i)<=chromEnd) {
												final char b1 = Character.toUpperCase(bases.charAt(readpos+i));
												if(refpos+i==C_arg_pos) {
													base_at_27_135_314 = b1;
													}
												sb.append(b1);
												}
											}
										refpos+=len;
										readpos+=len;

										break;
									default: throw new IllegalStateException(op.toString());
								}//end switch
							}//end loop cigar
							final String hisDNASeq = sb.toString();
							final String revcomp = AcidNucleics.reverseComplement(hisDNASeq);
							final StringBuilder pep = new StringBuilder();
							for(int i=0;i+2< revcomp.length();i+=3) {
								pep.append(geneticCode.translate(
									revcomp.charAt(i),
									revcomp.charAt(i+1),
									revcomp.charAt(i+2)));
								}
							final StringBuilder symbol = new StringBuilder();
							int i=0;
							while(i<pep.length()) {
								int n=0;
								char c= pep.charAt(i);
								while(i<pep.length() && pep.charAt(i)==c) {
									i++;
									n++;
									}
								final AminoAcids.AminoAcid aa = AminoAcids.getAminoAcidFromOneLetterCode(c);
								symbol.append(n).append(aa==null?"***":aa.getThreeLettersCode());
								}
							counter.incr(sample+"\t"+bamPath+"\t"+insertion_count+"\t"+deletion_count+"\t"+base_at_27_135_314+"\t"+hisDNASeq+"\t"+pep+"\t"+symbol);
							}//end while / iter
						}//end iterator
						int n=0;
						for(final String k: counter.keySetDecreasing()) {
							if(counter.count(k)<= this.min_occurence) continue;
							System.out.println(k+"\t"+counter.count(k)+"\t"+(counter.count(k)/(double)counter.getTotal()));
							n++;
							if(n>=2) break;
							}
						if(System.out.checkError()) {
							LOG.warn("ABORT");
							break;
							}
					}//end for SamReader
				}//end for path
			
			return 0;
		} catch(final Throwable err) {
			LOG.error(err);
			return -1;
		}
	}
	
	public static void main(String[] args) {
		new ScanHoxa1().instanceMainWithExit(args);
	}
	
}
