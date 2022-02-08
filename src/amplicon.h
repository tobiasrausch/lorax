#ifndef AMPLICON_H
#define AMPLICON_H

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/icl/split_interval_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "bed.h"
#include "var.h"
#include "util.h"

namespace lorax
{

  struct AmpliconConfig {
    uint16_t minQual;
    uint32_t minClip;
    uint32_t bpuncertain;
    std::string sample;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
    boost::filesystem::path bedfile;
    boost::filesystem::path vcffile;
  };

  template<typename TConfig>
  inline void
  outputReads(TConfig const& c, std::set<std::size_t> const& candidates) {
    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    // Parse BAM
    int32_t oldId = -1;
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) continue;
      if (oldId != rec->core.tid) {
	oldId = rec->core.tid;
	if ((oldId >= 0) && (oldId < hdr->n_targets)) {
      	  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Scanning " << hdr->target_name[oldId] << " for primary sequence" << std::endl;
	}
      }
      std::size_t seed = hash_string(bam_get_qname(rec));
      if (candidates.find(seed) == candidates.end()) continue;

      // Get read sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      
      // Output read
      dataOut << ">" << bam_get_qname(rec) << std::endl;
      dataOut << sequence << std::endl;
    }
    // Clean-up
    bam_destroy1(rec);
    
    // Close output file
    dataOut.pop();
    dataOut.pop();
    
    // Clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }
      
    
  template<typename TConfig>
  inline int32_t
  selectReads(TConfig const& c, std::set<std::size_t>& candidates) {
    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Load bcf file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);

    // Parse amplicons
    typedef boost::icl::interval_set<uint32_t> TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome scanRegions;
    int32_t nreg = _parseBedIntervals(c.bedfile.string(), hdr, scanRegions);
    if (nreg == 0) {
      std::cerr << "Warning: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
      return 1;
    } else {
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Parsed " << nreg << " amplicon regions." << std::endl;
    }

    // Discarded reads
    std::set<std::size_t> discarded;
    
    // Breakpoint support
    typedef std::set<uint32_t> TBpLocation;
    typedef std::map<std::size_t, TBpLocation> TReadClip;
    TReadClip readClip;
    
    // Assign reads to SNPs
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      std::string chrName(hdr->target_name[refIndex]);

      // Any amplicon regions?
      if (scanRegions[refIndex].empty()) continue;
      
      // Load het. markers
      typedef std::vector<BiallelicVariant> TPhasedVariants;
      TPhasedVariants pv;
      if (!_loadVariants(ibcffile, bcfidx, bcfhdr, c.sample, chrName, pv)) continue;

      // Sort variants
      std::sort(pv.begin(), pv.end(), SortVariants<BiallelicVariant>());

      // Load reference
      int32_t seqlen = -1;
      char* seq = NULL;
      seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);    

      // Assign reads to haplotypes
      for(typename TChrIntervals::iterator it = scanRegions[refIndex].begin(); it != scanRegions[refIndex].end(); ++it) {
	if ((it->lower() < it->upper()) && (it->lower() >= 0) && (it->upper() <= hdr->target_len[refIndex])) {
	  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing " << hdr->target_name[refIndex] << ':' << it->lower() << '-' << it->upper() << std::endl;

	  // Local haplotype assignment
	  std::set<std::size_t> h1;
	  std::set<std::size_t> h2;

	  // Breakpoint support
	  int32_t leftstart = 0;
	  int32_t leftend = hdr->target_len[refIndex];
	  if (c.bpuncertain < it->lower()) leftstart = it->lower() - c.bpuncertain;
	  if (it->lower() + c.bpuncertain < hdr->target_len[refIndex]) leftend = it->lower() + c.bpuncertain;
	  std::vector<uint16_t> leftBp((leftend - leftstart), 0);      
	  int32_t rightstart = 0;
	  int32_t rightend = hdr->target_len[refIndex];
	  if (c.bpuncertain < it->upper()) rightstart = it->upper() - c.bpuncertain;
	  if (it->upper() + c.bpuncertain < hdr->target_len[refIndex]) rightend = it->upper() + c.bpuncertain;
	  std::vector<uint16_t> rightBp((rightend - rightstart), 0);      
	  
	  // Scan region
	  hts_itr_t* iter = sam_itr_queryi(idx, refIndex, it->lower(), it->upper());
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile, iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minQual) || (rec->core.tid<0)) continue;
	    std::string qname = bam_get_qname(rec);
	    std::size_t seed = hash_string(bam_get_qname(rec));


	    // Split reads
	    {
	      uint32_t* cigar = bam_get_cigar(rec);
	      int32_t gp = rec->core.pos; // Genomic position
	      int32_t sp = 0; // Sequence position
	      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
		  if (bam_cigar_oplen(cigar[i]) > c.minClip) {
		    if (sp == 0) { // Left breakpoint
		      if (std::abs(gp - (int32_t) it->lower()) < c.bpuncertain) {
			if (readClip.find(seed) == readClip.end()) readClip[seed] = TBpLocation();
			readClip[seed].insert(it->lower());
			if ((gp >= leftstart) && (gp < leftend)) ++leftBp[gp-leftstart];
		      } else if (std::abs(gp - (int32_t) it->upper()) < c.bpuncertain) {
			if (readClip.find(seed) == readClip.end()) readClip[seed] = TBpLocation();
			readClip[seed].insert(it->lower());
			if ((gp >= rightstart) && (gp < rightend)) ++rightBp[gp-rightstart];
		      }
		    }
		  }
		  if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
		}
		else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
		else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
		else if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		  gp += bam_cigar_oplen(cigar[i]);
		  sp += bam_cigar_oplen(cigar[i]);
		} else {
		  std::cerr << "Unknown Cigar options" << std::endl;
		  return 1;
		}
	      }
	    }
	    

	    // Amplicon reads
	    bool ampliconRead = true;
	    if ((rec->core.pos < it->lower()) && (it->lower() - rec->core.pos > c.bpuncertain)) ampliconRead = false;
	    uint32_t alignEnd = lastAlignedPosition(rec);
	    if ((alignEnd > it->upper()) && (alignEnd - it->upper() > c.bpuncertain)) ampliconRead = false;
	    if (ampliconRead) {
	      uint32_t hp1votes = 0;
	      uint32_t hp2votes = 0;
	      int32_t regstart = rec->core.pos;
	      if (regstart < (int) it->lower()) regstart = it->lower();
	      TPhasedVariants::iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(regstart), SortVariants<BiallelicVariant>());
	      int32_t regend = alignEnd;
	      if (regend > (int) it->upper()) regend = it->upper();
	      TPhasedVariants::iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(regend), SortVariants<BiallelicVariant>());

	      // Any variants?
	      if (vIt != vItEnd) {
		// Get read sequence
		std::string sequence;
		sequence.resize(rec->core.l_qseq);
		uint8_t* seqptr = bam_get_seq(rec);
		for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	    
		// Parse CIGAR
		uint32_t* cigar = bam_get_cigar(rec);
		for(;vIt != vItEnd; ++vIt) {
		  int32_t gp = rec->core.pos; // Genomic position
		  int32_t sp = 0; // Sequence position
		  bool varFound = false;
		  for (std::size_t i = 0; ((i < rec->core.n_cigar) && (!varFound)); ++i) {
		    if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
		    else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {}
		    else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
		    else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
		    else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
		    else if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		      if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
			gp += bam_cigar_oplen(cigar[i]);
			sp += bam_cigar_oplen(cigar[i]);
		      } else {
			for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
			  if (gp == vIt->pos) {
			    varFound = true;
			    // Check REF allele
			    if (vIt->ref == seq[gp]) {
			      if (sp < (int) sequence.size()) {
				if ((sequence[sp] == vIt->alt) && (sequence[sp] != vIt->ref)) {
				  // ALT supporting read
				  ++vIt->asup;
				  if (vIt->hap) ++hp1votes;
				  else ++hp2votes;
				} else if ((sequence[sp] != vIt->alt) && (sequence[sp] == vIt->ref)) {
				  // REF supporting read
				  ++vIt->rsup;
				  if (vIt->hap) ++hp2votes;
				  else ++hp1votes;
				}
			      }
			    }
			  }
			}
		      }
		    } else {
		      std::cerr << "Unknown Cigar options" << std::endl;
		      return 1;
		    }
		  }
		}
	      }
	      int32_t hp = 0;
	      if (hp1votes > 2*hp2votes) hp = 1;
	      else if (hp2votes > 2*hp1votes) hp = 2;
	      if (hp) {
		if (hp == 1) h1.insert(seed);
		else if (hp == 2) h2.insert(seed);
	      }
	    }
	  }
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	  
	  // Identify haplotype support
	  int32_t h1sup = 0;
	  int32_t h2sup = 0;
	  TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(it->lower()), SortVariants<BiallelicVariant>());
	  TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(it->upper()), SortVariants<BiallelicVariant>());	  
	  for(;vIt != vItEnd; ++vIt) {
	    // Debug code
	    //std::cerr << hdr->target_name[refIndex] << ':' << it->lower() << '-' << it->upper() << ':' << vIt->pos << ',' << vIt->ref << ',' << vIt->alt << ',' << vIt->hap << ':' << vIt->rsup << ',' << vIt->asup << std::endl;
	    if (vIt->hap) {
	      h1sup += vIt->asup;
	      h2sup += vIt->rsup;
	    } else {
	      h1sup += vIt->rsup;
	      h2sup += vIt->asup;
	    }
	  }

	  // Select amplicon haplotype
	  if (h1sup + h2sup > 0) {
	    if (h1sup > 2 * h2sup) {
	      candidates.insert(h1.begin(), h1.end());
	      discarded.insert(h2.begin(), h2.end());
	    } else if (h2sup > 2 * h1sup) {
	      discarded.insert(h1.begin(), h1.end());
	      candidates.insert(h2.begin(), h2.end());
	    }
	  }
	  
	  // Refine breakpoints
	  int32_t bestLeft = it->lower() - leftstart;
	  for(uint32_t i = 0; i < leftBp.size(); ++i) {
	    if (leftBp[i] > leftBp[bestLeft]) bestLeft = i;
	  }
	  int32_t bestRight = it->upper() - rightstart;
	  for(uint32_t i = 0; i < rightBp.size(); ++i) {
	    if (rightBp[i] > rightBp[bestRight]) bestRight = i;
	  }
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Estimated  " << hdr->target_name[refIndex] << ':' << (bestLeft + leftstart) << '-' << (bestRight + rightstart) << std::endl;
	}
      }
      if (seq != NULL) free(seq);
    }
    fai_destroy(fai);

    // Collect split-reads connecting at least 2 amplicons
    uint32_t concordant = 0;
    uint32_t discordant = 0;
    uint32_t unclear = 0;
    for(TReadClip::iterator rcit = readClip.begin(); rcit != readClip.end(); ++rcit) {
      if (rcit->second.size() > 1) {
	if (candidates.find(rcit->first) != candidates.end()) ++concordant;
	else if (discarded.find(rcit->first) != discarded.end()) ++discordant;
	else {
	  ++unclear;
	  candidates.insert(rcit->first);
	}
      }
    }
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Concordant: " << concordant << ", discordant: " << discordant << ", unclear: " << unclear << std::endl;

      
    // Close BCF
    bcf_hdr_destroy(bcfhdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ibcffile);
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    return candidates.size();
  }
    

  template<typename TConfig>
  inline int32_t
  runAmplicon(TConfig const& c) {
    
#ifdef PROFILE
    ProfilerStart("lorax.prof");
#endif

    // Select amplicon reads
    std::set<std::size_t> candidates;
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Selecting amplicon reads." << std::endl;
    int32_t candsize = selectReads(c, candidates);
    
    // Searching primary alignments for full sequence
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Fetching sequence information for " << candsize << " amplicon reads." << std::endl;
    outputReads(c, candidates);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }


  
  int amplicon(int argc, char** argv) {
    AmpliconConfig c;
    
    // Parameter
    boost::program_options::options_description generic("Options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. sequence quality")
      ("minclip,c", boost::program_options::value<uint32_t>(&c.minClip)->default_value(100), "min. clipping length")
      ("uncertain,u", boost::program_options::value<uint32_t>(&c.bpuncertain)->default_value(1000), "breakpoint uncertainty (in bp)")
      ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name (as in VCF/BCF file)")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file")
      ("bedfile,b", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "amplicon regions in BED format")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.fa.gz"), "gzipped output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.tumor), "input file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("genome")) || (!vm.count("vcffile")) || (!vm.count("bedfile"))) {
      std::cout << "Usage: lorax " << argv[0] << " [OPTIONS] -g <ref.fa> -b amplicons.bed -s NA12878 -v <snps.bcf> <tumor.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Check input BAM file
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.tumor) && boost::filesystem::is_regular_file(c.tumor) && boost::filesystem::file_size(c.tumor))) {
	std::cerr << "Input BAM file is missing: " << c.tumor.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.tumor.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.tumor.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.tumor.string() << std::endl;
	return 1;
      }
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
  
    // Check VCF/BCF file
    if (vm.count("vcffile")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
	return 1;
      }
      hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
      if (bcfidx == NULL) {
	std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
      hts_idx_destroy(bcfidx);
      bcf_close(ifile);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "lorax ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runAmplicon(c);
  }

}

#endif
