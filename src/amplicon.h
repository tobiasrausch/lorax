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
    uint32_t wincov;
    uint32_t bpuncertain;
    uint32_t ploidy;
    std::string sample;
    std::string outprefix;
    boost::filesystem::path genome;
    boost::filesystem::path tumor;
    boost::filesystem::path bedfile;
    boost::filesystem::path vcffile;
  };


  template<typename TConfig, typename TScanRegions, typename TEdges>
  inline void
  outputAmplicons(TConfig const& c, TScanRegions const& scanRegions, std::vector<uint32_t> const& components, TEdges& es, uint64_t const avgcov) {
    // Compute node degree
    std::vector<uint32_t> degree(components.size(), 0);
    for(uint32_t i = 0; i < components.size(); ++i) {
      for(uint32_t j = i + 1; j < components.size(); ++j) {
	if (es.find(std::make_pair(i, j)) != es.end()) {
	  degree[i] += es[std::make_pair(i, j)];
	  degree[j] += es[std::make_pair(i, j)];
	}
      }
    }
    
    // Output segments
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    std::string filename = c.outprefix + ".bed";
    std::ofstream ofile(filename.c_str());
    ofile << "chr\tstart\tend\tnodeid\test_start\test_end\tncontent\tdegree\testcn\tclusterid\tedges" << std::endl;
    for(uint32_t refIndex = 0; refIndex < scanRegions.size(); ++refIndex) {
      for(uint32_t i = 0; i < scanRegions[refIndex].size(); ++i) {
	// Estimate coverage
	double estcov = (double) c.ploidy * ((((double) c.wincov / (double) (scanRegions[refIndex][i].end - scanRegions[refIndex][i].start)) * (double) scanRegions[refIndex][i].cov) / (double) avgcov);
	
	ofile << hdr->target_name[refIndex] << '\t' << scanRegions[refIndex][i].start << '\t' << scanRegions[refIndex][i].end << '\t';
	ofile << scanRegions[refIndex][i].cid << "[label=\"" << hdr->target_name[refIndex] << ':' << scanRegions[refIndex][i].start << '-' << scanRegions[refIndex][i].end << "(" << scanRegions[refIndex][i].cid << ")" <<  "\"];" << '\t';
	ofile << scanRegions[refIndex][i].estStart << '\t';
	ofile << scanRegions[refIndex][i].estEnd << '\t';
	ofile << scanRegions[refIndex][i].nsum << '\t';
	ofile << degree[scanRegions[refIndex][i].cid] << "\t" << estcov << "\t" << components[scanRegions[refIndex][i].cid] << '\t';
	for(uint32_t id2 = scanRegions[refIndex][i].cid; id2 < components.size(); ++id2) {
	  if (es.find(std::make_pair(scanRegions[refIndex][i].cid, id2)) != es.end()) {
	    ofile << scanRegions[refIndex][i].cid << "--" << id2 << "[label=\"" << es[std::make_pair(scanRegions[refIndex][i].cid, id2)] << "\"];";
	  }
	}
	ofile << std::endl;
      }
    }
    ofile.close();
    bam_hdr_destroy(hdr);
    sam_close(samfile);
  }
  
  template<typename TConfig>
  inline void
  outputReads(TConfig const& c, std::set<std::size_t> const& candidates) {
    std::string filename = c.outprefix + ".reads";
    std::ofstream ofile(filename.c_str());
    for(std::set<std::size_t>::const_iterator it = candidates.begin(); it != candidates.end(); ++it) {
      ofile << *it << std::endl;
    }
    ofile.close();
  }


  template<typename TEdgeSupport>
  inline void
  ampconnect(TEdgeSupport const& es, std::vector<uint32_t>& components) {
    for(typename TEdgeSupport::const_iterator ie = es.begin(); ie != es.end(); ++ie) {
      if (components[ie->first.first] != components[ie->first.second]) {
	uint32_t oldid = components[ie->first.second];
	for(uint32_t i = 0; i < components.size(); ++i) {
	  if (components[i] == oldid) components[i] = components[ie->first.first];
	}
      }
    }
  }
  
  template<typename TReadAmplicons, typename TEdgeSupport>
  inline void
  computeBridges(TReadAmplicons const& readAmp, TEdgeSupport& es) {
    typedef typename TReadAmplicons::mapped_type TAmpliconIds;
    for(typename TReadAmplicons::const_iterator rcit = readAmp.begin(); rcit != readAmp.end(); ++rcit) {
      for(typename TAmpliconIds::const_iterator it1 = rcit->second.begin(); it1 != rcit->second.end(); ++it1) {
	typename TAmpliconIds::const_iterator it2 = it1;
	++it2;
	for(; it2 != rcit->second.end(); ++it2) {
	  uint32_t id1 = *it1;
	  uint32_t id2 = *it2;
	  if (id1 > id2) {
	    uint32_t tmp = id1;
	    id1 = id2;
	    id2 = tmp;
	  }
	  if (es.find(std::make_pair(id1, id2)) == es.end()) es.insert(std::make_pair(std::make_pair(id1, id2), 1));
	  else ++es[std::make_pair(id1, id2)];
	}
      }
    }
  }
  
    
  template<typename TConfig, typename TScanRegions, typename TReadAmplicons>
  inline int32_t
  selectReads(TConfig const& c, TScanRegions& scanRegions, std::set<std::size_t>& candidates, TReadAmplicons& readAmp) {
    typedef typename TReadAmplicons::mapped_type TAmpliconIds;
    
    // Open file handles
    samFile* samfile = sam_open(c.tumor.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.tumor.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Load bcf file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);

    // Discarded reads
    std::set<std::size_t> discarded;
    
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
      std::sort(pv.begin(), pv.end());

      // Load reference
      int32_t seqlen = -1;
      char* seq = NULL;
      seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr->target_len[refIndex], 0);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      }

      // Assign reads to haplotypes
      for(uint32_t ri = 0; ri < scanRegions[refIndex].size(); ++ri) {
	if ((scanRegions[refIndex][ri].start < scanRegions[refIndex][ri].end) && (scanRegions[refIndex][ri].start >= 0) && (scanRegions[refIndex][ri].end <= hdr->target_len[refIndex])) {
	  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Processing " << hdr->target_name[refIndex] << ':' << scanRegions[refIndex][ri].start << '-' << scanRegions[refIndex][ri].end << std::endl;

	  // Set N-content
	  uint32_t nsum = 0;
	  for(uint32_t k = scanRegions[refIndex][ri].start; k < scanRegions[refIndex][ri].end; ++k) {
	    if (nrun[k]) ++nsum;
	  }
	  scanRegions[refIndex][ri].nsum = nsum;
	  
	  // Local haplotype assignment
	  std::set<std::size_t> h1;
	  std::set<std::size_t> h2;

	  // Breakpoint support
	  int32_t leftstart = 0;
	  int32_t leftend = hdr->target_len[refIndex];
	  if (c.bpuncertain < scanRegions[refIndex][ri].start) leftstart = scanRegions[refIndex][ri].start - c.bpuncertain;
	  if (scanRegions[refIndex][ri].start + c.bpuncertain < hdr->target_len[refIndex]) leftend = scanRegions[refIndex][ri].start + c.bpuncertain;
	  std::vector<uint16_t> leftBp((leftend - leftstart), 0);      
	  int32_t rightstart = 0;
	  int32_t rightend = hdr->target_len[refIndex];
	  if (c.bpuncertain < scanRegions[refIndex][ri].end) rightstart = scanRegions[refIndex][ri].end - c.bpuncertain;
	  if (scanRegions[refIndex][ri].end + c.bpuncertain < hdr->target_len[refIndex]) rightend = scanRegions[refIndex][ri].end + c.bpuncertain;
	  std::vector<uint16_t> rightBp((rightend - rightstart), 0);      
	  
	  // Scan region
	  hts_itr_t* iter = sam_itr_queryi(idx, refIndex, scanRegions[refIndex][ri].start, scanRegions[refIndex][ri].end);
	  bam1_t* rec = bam_init1();
	  uint64_t covsum = 0;
	  while (sam_itr_next(samfile, iter, rec) >= 0) {
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	    if ((rec->core.qual < c.minQual) || (rec->core.tid<0)) continue;
	    std::string qname = bam_get_qname(rec);
	    std::size_t seed = hash_string(bam_get_qname(rec));


	    // Split reads and coverage
	    {
	      uint32_t* cigar = bam_get_cigar(rec);
	      int32_t gp = rec->core.pos; // Genomic position
	      int32_t sp = 0; // Sequence position
	      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		if ((bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) || (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP)) {
		  if (bam_cigar_oplen(cigar[i]) > c.minClip) {
		    if (sp == 0) { // Left breakpoint
		      if (std::abs(gp - (int32_t) scanRegions[refIndex][ri].start) < c.bpuncertain) {
			if (readAmp.find(seed) == readAmp.end()) readAmp[seed] = TAmpliconIds();
			readAmp[seed].insert(scanRegions[refIndex][ri].cid);
			if ((gp >= leftstart) && (gp < leftend)) ++leftBp[gp-leftstart];
		      }
		    } else {
		      if (std::abs(gp - (int32_t) scanRegions[refIndex][ri].end) < c.bpuncertain) {
			if (readAmp.find(seed) == readAmp.end()) readAmp[seed] = TAmpliconIds();
			readAmp[seed].insert(scanRegions[refIndex][ri].cid);
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
		  for(uint32_t k = 0; k < bam_cigar_oplen(cigar[i]); ++k, ++gp, ++sp) {
		    if ((gp >= (int) scanRegions[refIndex][ri].start) && (gp < (int) scanRegions[refIndex][ri].end)) ++covsum;
		  }
		} else {
		  std::cerr << "Unknown Cigar options" << std::endl;
		  return 1;
		}
	      }
	    }
	    

	    // Amplicon reads
	    bool ampliconRead = true;
	    if ((rec->core.pos < scanRegions[refIndex][ri].start) && (scanRegions[refIndex][ri].start - rec->core.pos > c.bpuncertain)) ampliconRead = false;
	    uint32_t alignEnd = lastAlignedPosition(rec);
	    if ((alignEnd > scanRegions[refIndex][ri].end) && (alignEnd - scanRegions[refIndex][ri].end > c.bpuncertain)) ampliconRead = false;
	    if (ampliconRead) {
	      uint32_t hp1votes = 0;
	      uint32_t hp2votes = 0;
	      int32_t regstart = rec->core.pos;
	      if (regstart < (int) scanRegions[refIndex][ri].start) regstart = scanRegions[refIndex][ri].start;
	      TPhasedVariants::iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(regstart));
	      int32_t regend = alignEnd;
	      if (regend > (int) scanRegions[refIndex][ri].end) regend = scanRegions[refIndex][ri].end;
	      TPhasedVariants::iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(regend));

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
	  scanRegions[refIndex][ri].cov = covsum;
	  
	  // Identify haplotype support
	  int32_t h1sup = 0;
	  int32_t h2sup = 0;
	  TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(scanRegions[refIndex][ri].start));
	  TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(scanRegions[refIndex][ri].end));
	  for(;vIt != vItEnd; ++vIt) {
	    // Debug code
	    //std::cerr << hdr->target_name[refIndex] << ':' << scanRegions[refIndex][ri].start << '-' << scanRegions[refIndex][ri].end << ':' << vIt->pos << ',' << vIt->ref << ',' << vIt->alt << ',' << vIt->hap << ':' << vIt->rsup << ',' << vIt->asup << std::endl;
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
	  int32_t bestLeft = scanRegions[refIndex][ri].start - leftstart;
	  for(uint32_t i = 0; i < leftBp.size(); ++i) {
	    if (leftBp[i] > leftBp[bestLeft]) bestLeft = i;
	  }
	  int32_t bestRight = scanRegions[refIndex][ri].end - rightstart;
	  for(uint32_t i = 0; i < rightBp.size(); ++i) {
	    if (rightBp[i] > rightBp[bestRight]) bestRight = i;
	  }
	  scanRegions[refIndex][ri].estStart = bestLeft + leftstart;
	  scanRegions[refIndex][ri].estEnd = bestRight + rightstart;
	}
      }
      if (seq != NULL) free(seq);
    }
    fai_destroy(fai);

    // Collect split-reads connecting at least 2 amplicons
    uint32_t concordant = 0;
    uint32_t discordant = 0;
    uint32_t unclear = 0;
    for(typename TReadAmplicons::iterator rcit = readAmp.begin(); rcit != readAmp.end(); ++rcit) {
      if (rcit->second.size() > 1) {
	// True split
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

    // Parse amplicons
    typedef std::vector<BedEntry> TChrIntervals;
    typedef std::vector<TChrIntervals> TRegionsGenome;
    TRegionsGenome scanRegions;
    int32_t nreg = _parseBedIntervals(c, scanRegions);
    if (nreg == 0) {
      std::cerr << "Error: Couldn't parse BED intervals. Do the chromosome names match?" << std::endl;
      return 1;
    } else {
      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Parsed " << nreg << " amplicon regions." << std::endl;
    }

    // Estimate coverage
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Estimate coverage." << std::endl;
    uint64_t avgcov = 0;
    uint64_t sdcov = 0;
    estimateCoverage(c.tumor.string(), c.genome.string(), c.wincov, avgcov, sdcov);
    if (avgcov == 0) {
      std::cerr << "Error: Zero background coverage!" << std::endl;
      return 1;
    }
    
    // Select amplicon reads
    std::set<std::size_t> candidates;
    typedef std::set<uint32_t> TAmpliconIds;
    typedef std::map<std::size_t, TAmpliconIds> TReadAmplicons;
    TReadAmplicons readAmp;
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Selecting amplicon reads." << std::endl;
    int32_t candsize = selectReads(c, scanRegions, candidates, readAmp);

    // Amplicon bridges
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Computing amplicon bridges" << std::endl;
    typedef std::pair<uint32_t, uint32_t> TEdge;
    typedef std::map<TEdge, uint32_t> TEdgeSupport;
    TEdgeSupport es;
    computeBridges(readAmp, es);

    // Connected components (should be one)
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Computing connected components" << std::endl;
    std::vector<uint32_t> components(nreg);
    for(uint32_t i = 0; i < components.size(); ++i) components[i] = i;
    ampconnect(es, components);

    // Output amplicons
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Output amplicon connections." << std::endl;
    outputAmplicons(c, scanRegions, components, es, avgcov);
      
    // Searching primary alignments for full sequence
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Output " << candsize << " amplicon reads." << std::endl;
    outputReads(c, candidates);
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    // End
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Done." << std::endl;
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
      ("wincov,n", boost::program_options::value<uint32_t>(&c.wincov)->default_value(1000), "coverage window length")
      ("uncertain,u", boost::program_options::value<uint32_t>(&c.bpuncertain)->default_value(1000), "breakpoint uncertainty (in bp)")
      ("ploidy,p", boost::program_options::value<uint32_t>(&c.ploidy)->default_value(2), "ploidy")
      ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name (as in VCF/BCF file)")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input VCF/BCF file")
      ("bedfile,b", boost::program_options::value<boost::filesystem::path>(&c.bedfile), "amplicon regions in BED format")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome fasta file")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
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
