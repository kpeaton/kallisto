#include "enhancedoutput.h"

HighResTimer::HighResTimer()
{
#if defined(_MSC_VER)
	QueryPerformanceFrequency(&frequency);
#endif
	resetTimer();
}

void HighResTimer::resetTimer()
{
#if defined(_MSC_VER)
	QueryPerformanceCounter(&reset_time);
#else
	reset_time = std::chrono::high_resolution_clock::now();
#endif
	previous_time = reset_time;
}

HighResTimer::duration HighResTimer::timeSinceReset()
{
#if defined(_MSC_VER)
	QueryPerformanceCounter(&current_time);
	duration elapsed = std::chrono::microseconds(((current_time.QuadPart - reset_time.QuadPart) * 1000000) / frequency.QuadPart);
#else
	current_time = std::chrono::high_resolution_clock::now();
	duration elapsed = current_time - reset_time;
#endif
	previous_time = current_time;
	return elapsed;
}

HighResTimer::duration HighResTimer::timeSincePrevious()
{
#if defined(_MSC_VER)
	QueryPerformanceCounter(&current_time);
	duration elapsed = std::chrono::microseconds(((current_time.QuadPart - previous_time.QuadPart) * 1000000) / frequency.QuadPart);
#else
	current_time = std::chrono::high_resolution_clock::now();
	duration elapsed = current_time - previous_time;
#endif
	previous_time = current_time;
	return elapsed;
}

EnhancedOutput::EnhancedOutput(KmerIndex &index, const ProgramOptions& opt)
	: enhancedoutput(opt.pseudobam && !opt.exon_coords_file.empty()),
	  sortedbam(opt.sortedbam),
	  outputunmapped(false),
	  outputbed(opt.pseudobam && !opt.exon_coords_file.empty() && !opt.bed_file.empty())
{
	if (enhancedoutput) {

		// Load the exon/intron coordinate map data
		std::ifstream file(opt.exon_coords_file);
		CSVRow row;
		std::string last_key;
		while (file >> row) {
			if (last_key == row[0]){
				std::get<6>(exon_map[last_key]).emplace_back(std::vector<int>({ std::stoi(row[2]), std::stoi(row[3]), std::stoi(row[4]) }));
			} else {
				IntronFlag intron_flag = intronNone;
				int segment_start = std::stoi(row[2]);
				int segment_end = std::stoi(row[3]);
				int genome_position = std::stoi(row[4]);
				int pair_start = -1;
				int pair_end = -1;
				if (row[0].back() == ')') {
					if ((last_key.back() == ')') && (last_key.substr(0, last_key.find("::") + 2) == row[0].substr(0, row[0].find("::") + 2))) {  // Store pairs together?!
						intron_flag = intronEnd;
						int start_coord = std::get<6>(exon_map[last_key])[0][2];
						int end_coord = genome_position + segment_end - segment_start;
						pair_start = start_coord + 39;
						pair_end = start_coord + 59;
						std::get<3>(exon_map[last_key]) = intronStart;
						std::get<4>(exon_map[last_key]) = end_coord - 60;
						std::get<5>(exon_map[last_key]) = end_coord - 40;
					} else {
						intron_flag = intronFull;
					}
				}
				last_key = row[0];
				exon_map.emplace(last_key, std::make_tuple(row[5], row[6], std::stoi(row[1]), intron_flag, pair_start, pair_end, std::vector<std::vector<int>>({ std::vector<int>({ segment_start, segment_end, genome_position }) })));
			}
		}
		file.close();

		// Collect information for the bam header
		for (int i = 0; i < index.num_trans; i++) {

			const char * key = index.target_names_[i].c_str();
			ExonMap::const_iterator map_entry = exon_map.find(key);
			if (map_entry == exon_map.end()) {
				std::cerr << "Transcript name could not be found in exon coordinate file: " << key << std::endl;
				exit(1);
			}

			int strand = std::get<2>(map_entry->second);
			auto &entry = std::get<6>(map_entry->second);
			int len;
			if (strand < 0) {
				len = entry.front()[2] + entry.front()[1] - entry.front()[0];
			} else {
				len = entry.back()[2] + entry.back()[1] - entry.back()[0];
			}

			key = std::get<1>(map_entry->second).c_str();
			if (ref_seq_map.find(key) == ref_seq_map.end()) {
				ref_seq_map.emplace(key, std::vector<int>({ len, -1, 0 }));
			} else {
				ref_seq_map[key][0] = (std::max)(ref_seq_map[key][0], len);
			}

		}

		if (sortedbam) {

			// Build header and initialize output file streams for sorting
			std::ostringstream header_stream;
			header_stream << "@HD\tVN:1.0\tSO:coordinate\n";

			sort_dir = opt.output + "/sorting";
			mkdir(sort_dir.c_str(), 0777);  // Move this to CheckOptions?

			int refid = 0;
			std::fstream::openmode stream_flags = std::fstream::in | std::fstream::out | std::fstream::trunc | std::fstream::binary;
			for (auto &entry : ref_seq_map) {

				header_stream << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second[0] << "\n";
				sort_file_map[entry.first] = new std::fstream(sort_dir + "/" + entry.first, stream_flags); // Use std::move instead of new? Use emplace?
				entry.second[1] = refid++;

			}

			header_stream << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
			sam_header = header_stream.str();

		} else {

			// Write the bam header
			std::cout << "@HD\tVN:1.0\n";
			for (auto const &entry : ref_seq_map) {
				std::cout << "@SQ\tSN:" << entry.first << "\tLN:" << entry.second[0] << "\n";
			}
			std::cout << "@PG\tID:kallisto\tPN:kallisto\tVN:" << KALLISTO_VERSION << "\n";
			std::cout.flush();

		}

		if (outputbed) {
			bed_file = opt.bed_file;
		}

	} else if (opt.pseudobam) {

		// Write the bam header
		index.writePseudoBamHeader(std::cout);

	}
}

EnhancedOutput::EnhancedOutput()
	: enhancedoutput(false),
	  sortedbam(false),
	  outputunmapped(false),
	  outputbed(false)
{
}

void EnhancedOutput::processAlignment(std::string &ref_name, int& strand, int &posread, int &posmate, int slen1, int slen2, char *cig, std::vector<uint> &bam_cigar, uint &align_len)
{
	std::string trans_name = ref_name;
	ExonMap::const_iterator map_entry = exon_map.find(ref_name);
	ref_name = std::get<1>(map_entry->second);
	strand = (std::get<2>(map_entry->second) < 0) ? -1 : 1;  // Use trsense from findPosition instead?!
	bool negstrand = strand < 0;
	IntronFlag intron_flag = std::get<3>(map_entry->second);

	std::string cig_string;
	uint op_len;
	int read_rem = slen1;
	int mate_rem = (posmate == 0) ? 0 : slen2;
	int read_offset = 0;
	int start_coord = 0;
	int end_coord = 0;

	for (auto &span : std::get<6>(map_entry->second)) {

		if (read_rem > 0) {  // Not done mapping read

			if (read_rem < slen1) {  // In the process of mapping read

				// Find span of skipped intron:
				if (negstrand) {
					start_coord = span[2] + span[1] - span[0];
					end_coord = read_offset;
				} else {
					start_coord = read_offset;
					end_coord = span[2];
				}

				// Update CIGAR string:
				if (sortedbam) {
					buildBAMCigar(bam_cigar, align_len, negstrand, end_coord - start_coord - 1, 3);
				} else{
					buildSAMCigar(cig_string, negstrand, end_coord - start_coord - 1, 'N');
				}

				// Store BED file information:
				if (outputbed) {
					mapJunction(ref_name, trans_name, negstrand, start_coord, end_coord, 0, 0, -1, -1);
				}

				// Find overlap with exon segment:
				read_offset = span[1] - span[0] + 1;
				if (read_rem > read_offset){
					op_len = read_offset;
					read_offset = (negstrand) ? span[2] : read_offset + span[2] - 1;
				} else {
					op_len = read_rem;
					if (negstrand) {
						posread = start_coord - read_rem + 1;
					}
				}
				read_rem -= op_len;

				// Update CIGAR string:
				if (sortedbam) {
					buildBAMCigar(bam_cigar, align_len, negstrand, op_len, 0);
				} else{
					buildSAMCigar(cig_string, negstrand, op_len, 'M');
				}

			} else if (posread <= span[1]) {  // Begin mapping read

				// Find soft clipping at beginning of read:
				if (posread < span[0]) {
					op_len = span[0] - posread;
					read_rem -= op_len;
					if (sortedbam) {
						buildBAMCigar(bam_cigar, align_len, false, op_len, 4);
					} else{
						buildSAMCigar(cig_string, false, op_len, 'S');
					}
				}

				// Find overlap with exon segment:
				read_offset = posread + slen1 - span[1] - 1;
				if (read_offset > 0) {
					op_len = read_rem - read_offset;
					if (negstrand) {
						read_offset = span[2];
					} else {
						read_offset = span[2] + span[1] - span[0];
						posread += span[2] - span[0];
					}
				} else {
					op_len = read_rem;
					posread = (negstrand) ? span[2] - read_offset : posread + span[2] - span[0];
				}
				read_rem -= op_len;

				// Update CIGAR string:
				if (sortedbam) {
					buildBAMCigar(bam_cigar, align_len, negstrand, op_len, 0);
				} else{
					buildSAMCigar(cig_string, negstrand, op_len, 'M');
				}

			}

		}

		if (mate_rem > 0) {  // Not done mapping mate

			if (mate_rem < slen2) {  // In the process of mapping mate

				mate_rem -= span[1] - span[0] + 1;
				posmate = span[2] - mate_rem;

			} else if (posmate <= span[1]) {  // Begin mapping mate

				if (negstrand) {
					mate_rem = posmate + slen2 - span[1] - 1;
					posmate = span[2] - mate_rem;
				} else {
					mate_rem = 0;
					posmate += span[2] - span[0];
				}

			}

		}

		if ((read_rem <= 0) && (mate_rem <= 0)) {  // Both have been mapped
			break;
		}

	}

	if (read_rem > 0) {  // Account for read overhangs
		if (negstrand) {
			posread = read_offset - read_rem;
		}
		if (sortedbam) {
			buildBAMCigar(bam_cigar, align_len, negstrand, read_rem, 4);
		} else{
			buildSAMCigar(cig_string, negstrand, read_rem, 'S');
		}
	}

	if (mate_rem == slen2) {  // Account for mate completely outside segment (is this even necessary?!)
		std::cerr << "mate outside segment: " << trans_name << std::endl;
		auto &span = std::get<6>(map_entry->second).back();
		posmate = (negstrand) ? span[2] - posmate - slen2 + span[1] + 1 : posmate + span[2] - span[0];
	}

	if (outputbed && (intron_flag != intronNone)) {

		auto &span = std::get<6>(map_entry->second)[0];
		start_coord = span[2];
		end_coord = span[2] + span[1] - span[0];
		trans_name = trans_name.substr(0, trans_name.find("::")) + '-';

		switch (intron_flag) {  // Could perhaps make this more compact and less redundant
			case intronStart: {
				if ((posread >= start_coord) && (posread < start_coord + 50) &&
					(posread + slen1 >= start_coord + 50) && (posread + slen1 < end_coord) &&
					(posmate < end_coord)) {
					mapJunction(ref_name, trans_name + std::to_string(start_coord + 50), negstrand, start_coord + 39, start_coord + 59, 10, 10, std::get<4>(map_entry->second), std::get<5>(map_entry->second));
				}
				break;
			}
			case intronEnd: {
				if ((posread >= start_coord) && (posread < end_coord - 50) &&
					(posread + slen1 >= end_coord - 50) && (posread + slen1 < end_coord) &&
					(posmate + slen2 >= start_coord)) {
					mapJunction(ref_name, trans_name + std::to_string(end_coord - 50), negstrand, end_coord - 60, end_coord - 40, 10, 10, std::get<4>(map_entry->second), std::get<5>(map_entry->second));
				}
				break;
			}
			case intronFull: {
				if ((posread >= start_coord) && (posread < start_coord + 50) &&
					(posread + slen1 >= start_coord + 50) && (posread + slen1 < end_coord - 50) &&
					(posmate < end_coord)) {
					mapJunction(ref_name, trans_name + std::to_string(start_coord + 50), negstrand, start_coord + 39, start_coord + 59, 10, 10, end_coord - 60, end_coord - 40);
				}
				if ((posread >= start_coord + 50) && (posread < end_coord - 50) &&
					(posread + slen1 >= end_coord - 50) && (posread + slen1 < end_coord) &&
					(posmate + slen2 >= start_coord)) {
					mapJunction(ref_name, trans_name + std::to_string(end_coord - 50), negstrand, end_coord - 60, end_coord - 40, 10, 10, start_coord + 39, start_coord + 59);
				}
				break;
			}
		}
	}

	if (!sortedbam) {
		sprintf(cig, "%s", cig_string.c_str());
	}
}

void EnhancedOutput::buildSAMCigar(std::string &cig_string, bool prepend, uint op_len, const char cig_char)
{
	static char cig_[10];
	char *cig = &cig_[0];

	sprintf(cig, "%d%c", op_len, cig_char);  // Would to_string be faster?
	if (prepend) {
		cig_string.insert(0, cig);
	}
	else {
		cig_string += cig;
	}
}

void EnhancedOutput::buildBAMCigar(std::vector<uint> &bam_cigar, uint &align_len, bool prepend, uint op_len, uint cig_int)
{
	if (prepend) {
		bam_cigar.insert(bam_cigar.begin(), ((op_len << 4) | cig_int));
	} else {
		bam_cigar.push_back(((op_len << 4) | cig_int));
	}
	if (cig_int != 4) {
		align_len += op_len;
	}
}

void EnhancedOutput::mapJunction(std::string chrom_name, std::string trans_name, bool negstrand, int start_coord, int end_coord, int size1, int size2, int pair_start, int pair_end)
{
	JunctionKey key = std::make_tuple(chrom_name, start_coord, end_coord);
	if (junction_map.find(key) == junction_map.end()) {
		junction_map.emplace(key, std::make_tuple(trans_name, 1, (negstrand) ? '-' : '+', size1, size2, pair_start, pair_end));
	} else {
		std::get<1>(junction_map[key])++;
	}
}

void EnhancedOutput::outputBamAlignment(std::string ref_name, int posread, int flag, int slen, int posmate, int tlen, const char *n1, std::vector<uint> bam_cigar, uint align_len, const char *seq, const char *qual, int nmap, int strand, int id)
{
	uint *buffer = (uint*)(outBamBuffer);
	uint n_bytes = 0;
	uint name_len;
	uint mapq = 255;
	uint n_cigar;
	uint cigar_bytes;

	// refID
	buffer[1] = ref_seq_map[ref_name][1];

	// pos
	buffer[2] = posread - 1;

	// bin_mq_nl
	name_len = strlen(n1) + 1;
	buffer[3] = ((reg2bin(posread - 1, posread + align_len - 1)<<16) | (mapq<<8) | name_len);

	// flag_nc
	n_cigar = bam_cigar.size();
	buffer[4] = ((flag << 16) | n_cigar);

	// l_seq
	buffer[5] = slen;

	// next_refID
	buffer[6] = ref_seq_map[ref_name][1];

	// next_pos
	buffer[7] = posmate - 1;

	// tlen
	buffer[8] = tlen;

	n_bytes = 9 * sizeof(uint);

	// read_name
	memcpy(outBamBuffer + n_bytes, n1, name_len);
	n_bytes += name_len;

	// cigar
	cigar_bytes = n_cigar*sizeof(uint);
	memcpy(outBamBuffer + n_bytes, bam_cigar.data(), cigar_bytes);
	n_bytes += cigar_bytes;

	// seq
	packseq(seq, outBamBuffer + n_bytes, slen);
	n_bytes += (slen + 1)/2;

	// qual
	for (uint i = 0; i < slen; i++) {
		(outBamBuffer + n_bytes)[i] = qual[i] - 33;
	};
	n_bytes += slen;

	// attributes
	memcpy(outBamBuffer + n_bytes, "NHi", 3);
	memcpy(outBamBuffer + n_bytes + 3, &nmap, sizeof(int));
	n_bytes += 3 + sizeof(int);
	memcpy(outBamBuffer + n_bytes, (strand < 0) ? "XSA-": "XSA+", 4);      // Which to use?!
	//memcpy(outBamBuffer + n_bytes, bool(flag & 0x10) ? "XSA-" : "XSA+", 4);
	n_bytes += 4;

	// block_size
	buffer[0] = n_bytes - sizeof(uint);

	// Output to sorting file
	sort_file_map[ref_name]->write(outBamBuffer, n_bytes);
	ref_seq_map[ref_name][2]++; // Combine with sort_file_map?
}

void EnhancedOutput::outputSortedBam()
{
	std::fstream *sort_file;
	char *align_buffer;
	uint *sort_buffer;

	// Connect BAM output to stdout
	BGZF *bam_stream;
#ifdef _WIN32
	int result = _setmode(_fileno(stdout), _O_BINARY);
#endif
//	std::setvbuf(stdout, NULL, _IOFBF, 65536);
	bam_stream = bgzf_dopen(fileno(stdout), "w1");

	//// ABORT!!!
	//for (auto &entry : ref_seq_map) {
	//	sort_file_map[entry.first]->close();
	//	delete sort_file_map[entry.first];
	//	remove((sort_dir + "/" + entry.first).c_str());
	//}
	//bgzf_flush(bam_stream);
	//bgzf_close(bam_stream);
	//return;

	// Output header
	bgzf_write(bam_stream, "BAM\001", 4);
	int hlen = (int)sam_header.size();
	bgzf_write(bam_stream, (char*)&hlen, sizeof(hlen));
	bgzf_write(bam_stream, sam_header.c_str(), hlen);
	int nchr = (int)ref_seq_map.size();
	bgzf_write(bam_stream, (char*)&nchr, sizeof(nchr));
	for (auto const &entry : ref_seq_map) {
		int namelen = (int)(entry.first.size() + 1);
		int seqlen = entry.second[0];
		//std::cerr << namelen << ": " << entry.first.data() << " " << seqlen << std::endl;
		bgzf_write(bam_stream, (char*)&namelen, sizeof(namelen));
		bgzf_write(bam_stream, entry.first.data(), namelen);
		bgzf_write(bam_stream, (char*)&seqlen, sizeof(seqlen));
	}
	bgzf_flush(bam_stream);

	// Create buffers to store each set of alignments and sorting data
	uint max_file_size = 0;
	uint max_n_align = 0;
	for (auto &entry : ref_seq_map) {
		max_file_size = (std::max)(max_file_size, (uint)sort_file_map[entry.first]->tellg());
		max_n_align = (std::max)(max_n_align, (uint)entry.second[2]);
	}
	align_buffer = new char[max_file_size+1];
	sort_buffer = new uint[max_n_align * 2];

	// Sort alignments for each sorting file
	for (auto &entry : ref_seq_map) {

		// Get file stream
		sort_file = sort_file_map[entry.first];
		uint n_align = entry.second[2];
		//std::cerr << entry.first << ": " << n_align << std::endl;
		if (n_align == 0) {
			sort_file->close();
			delete sort_file;
			remove((sort_dir + "/" + entry.first).c_str());
			continue;
		}

		// Load data from file
		uint n_bytes = sort_file->tellg();  // uint64_t?? size_t??
		sort_file->seekg(std::fstream::beg);
		sort_file->read(align_buffer, n_bytes);
		sort_file->close();
		delete sort_file;
		remove((sort_dir + "/" + entry.first).c_str());

		// Collect sorting data
		for (uint position = 0, i = 0; i < n_align; i++) {
			uint *buffer = (uint*)(align_buffer + position);
			sort_buffer[i * 2] = buffer[2];
			sort_buffer[i * 2 + 1] = position;
			position += buffer[0] + sizeof(uint);
		}

		// Sort by genomic position
		qsort((void*) sort_buffer, n_align, sizeof(uint) * 2, funCompareArrays<uint, 2>);

		// Output sorted alignments
		for (uint i = 0; i < n_align; i++) {
			char *buffer = align_buffer + sort_buffer[i * 2 + 1];
			bgzf_write(bam_stream, buffer, *((uint*)buffer) + sizeof(uint));
		}

	}

	// Flush and close output and delete buffers
	bgzf_flush(bam_stream);
	bgzf_close(bam_stream);
	delete[] align_buffer;
	delete[] sort_buffer;
}

// calculate bin given an alignment covering [beg,end) (zero-based, half-close-half-open)
int EnhancedOutput::reg2bin(int beg, int end)
{
	--end;
	if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
	if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
	if (beg>>20 == end>>20) return ((1<<9)-1)/7  + (beg>>20);
	if (beg>>23 == end>>23) return ((1<<6)-1)/7  + (beg>>23);
	if (beg>>26 == end>>26) return ((1<<3)-1)/7  + (beg>>26);
	return 0;
}

// pack nucleotides for BAM
void EnhancedOutput::packseq(const char *seq_in, char *seq_out, uint seq_len)
{
	for (uint i = 0; i < seq_len/2; i++) {
		seq_out[i] = ((encodeNucleotide(seq_in[2 * i])<<4) | encodeNucleotide(seq_in[2 * i + 1]));
	};
	if (seq_len%2 == 1) {
		seq_out[seq_len/2] = (encodeNucleotide(seq_in[seq_len - 1])<<4);
	};
}

// encode nucleotides for packing
char EnhancedOutput::encodeNucleotide(char cc) // Use a std::map instead?
{
	switch (cc) {
		case ('=') : cc = 0; break;
		case ('A') : case ('a') : cc = 1; break;
		case ('C') : case ('c') : cc = 2; break;
		case ('M') : case ('m') : cc = 3; break;
		case ('G') : case ('g') : cc = 4; break;
		case ('R') : case ('r') : cc = 5; break;
		case ('S') : case ('s') : cc = 6; break;
		case ('V') : case ('v') : cc = 7; break;
		case ('T') : case ('t') : cc = 8; break;
		case ('W') : case ('w') : cc = 9; break;
		case ('Y') : case ('y') : cc = 10; break;
		case ('H') : case ('h') : cc = 11; break;
		case ('K') : case ('k') : cc = 12; break;
		case ('D') : case ('d') : cc = 13; break;
		case ('B') : case ('b') : cc = 14; break;
		case ('N') : case ('n') : cc = 15; break;
		default: cc = 15;
	};
	return cc;
};

std::string const& CSVRow::operator[](std::size_t index) const
{
	return m_data[index];
}

std::size_t CSVRow::size() const
{
	return m_data.size();
}

void CSVRow::readNextRow(std::istream& str)
{
	std::string line;
	std::getline(str, line);

	std::stringstream lineStream(line);
	std::string cell;

	m_data.clear();
	while (std::getline(lineStream, cell, ',')) {
		m_data.push_back(cell);
	}
	// This checks for a trailing comma with no data after it.
	if (!lineStream && cell.empty()) {
		// If there was a trailing comma then add an empty element.
		m_data.push_back("");
	}
}

std::istream& operator>>(std::istream& str, CSVRow& data)
{
	data.readNextRow(str);
	return str;
}
