#!/usr/bin/perl

use strict;

$| = 1;

use Net::FTP;

# Figure out whether to obtain a new copy of the RefSeq genomic database,
# use an existing one, or refresh existing data with new data from RefSeq.

my $newCopy = 0;

my $logDir = '.logs/setup';

system("mkdir -p $logDir");

&checkDataSource();

if ( $newCopy == 0 or $newCopy == 1 ) {
   
   # USE EXISTING genome database ($newCopy == 0).  This is essentially
   # a cleanup operation.
   # 
   # --> OR <--
   # 
   # OVERWRITE genome database with new data ($newCopy == 1): We've been
   # told to completely overwrite existing genome data with current RefSeq
   # data, treating this as a first-time installation.
   
   # 0. At this point, the data's already been downloaded and moved to the proper place.

   # 1. Rename directories as necessary.

   &renameNewRefSeqDirs('.genomeData');
   
   # 2. Generate metafiles describing the RefSeq data, using NCBI's taxonomy data to help.

   &generateRefSeqMetafiles();

   # 3. Extract FASTA sequences from the GenBank files as necessary.

   &genbankToFasta();
   
   # 4. Create training files for the ICMs as necessary.

   &createTrainingFiles();
   
   # 5. Build the local BLAST database from the RefSeq files.

   &generateBlastDB();

   # 6. Build/update the ICM suite.

   &buildICMs();
   
} elsif ( $newCopy == 2 ) {
   
   # REFRESH existing genome database with changes from RefSeq since last update.
   
   # We're only doing two things here: adding newly-deposited genomes,
   # and updating existing ones.  We're /not/ deleting genomes which aren't
   # in the current RefSeq instance (for several reasons, including preserving
   # manually-added organisms).  If you want old organisms removed, just delete
   # the subdirectory of .genomeData corresponding to the offending organism(s).
   
   # 0. At this point, the new RefSeq data has been extracted
   #    to '.tempGenomeData', and the existing genome data remains
   #    untouched in '.genomeData'.
   
   # 1 Rename genome directories - in the existing genome database - as necessary.

   &renameNewRefSeqDirs('.genomeData');

   # 2 Rename genome directories - in the temporary holding bin containing the new RefSeq data - as necessary.

   &renameNewRefSeqDirs('.tempGenomeData');
   
   # 3. Perform genome data updates as needed and delete the temp directory.
   
   &updateGenomeData();
   
   # 4. Generate metafiles describing the RefSeq data, using NCBI's taxonomy data to help.

   &generateRefSeqMetafiles();

   # 5. Extract FASTA sequences from the GenBank files as necessary.

   &genbankToFasta();
   
   # 6. Create training files for the ICMs as necessary.

   &createTrainingFiles();

   # 7. Build the local BLAST database from the RefSeq files.

   &generateBlastDB();
   
   # 8. Build/update the ICM suite.

   &buildICMs();
}








###########################################################
# SUBROUTINES
###########################################################

sub checkDataSource {
   
   if ( &checkRefSeqData() == 1 ) {
      
      print "\nLocal genome database detected.  Please choose one of the following:\n\n";

      print "Use (E)xisting genome data without changing it\n\n";

      print "(R)efresh database with newest RefSeq data (will not alter genomes which haven't been\n";
      print "changed in the RefSeq database since Phymm was installed)\n\n";

      print "completely (O)verwrite existing genome data with fresh/current copy from RefSeq\n\n";

      print "[Default is (E)]: ";

      my $response = <STDIN>;

      chomp $response;

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[eErRoO]/ ) ) {
	 
	 print "\nPlease enter E (default: use existing genome data), R (refresh with newest data) or O (overwrite with new copy): ";

	 $response = <STDIN>;

	 chomp $response;
      }
      
      print "\n";
      
      if ( $response =~ /[oO]/ ) {
	 
	 # Grab the newest data from the server.
	 
	 &downloadNewRefSeqData();
	 
	 # Wipe .genomeData.  If the .userAdded subdirectory's in there, this syntax should leave it alone.

	 system("rm -rf .genomeData/*");

	 # Extract downloaded tarball into .genomeData.
	 
	 &expandNewRefSeqData_genomeData();

	 $newCopy = 1;
	 
      } elsif ( $response =~ /[rR]/ ) {

	 # Refresh data with genomes added since install, and update genomes which've been altered since install.

	 &downloadNewRefSeqData();
	 
	 # Extract downloaded tarball into a temp directory for comparison/update.

	 &expandNewRefSeqData_tempDir();

	 $newCopy = 2;

      } else {
	 
	 # We're going to use existing data.

	 $newCopy = 0; # (reaffirming default value)
      }

   } else {
      
      # In this case, $checkRefSeqData() has returned 0 and there's no tarball and no data in .genomeData.  Do a fresh install.

      print "\nNo RefSeq genomic data detected; download now? [Y/n] ";

      my $response = <STDIN>;

      chomp $response;

      while ( length($response) > 1 or ( length($response) == 1 and $response !~ /[yYnN]/ ) ) {
	 
	 print "\nPlease enter Y (default: get a new copy of RefSeq genome data) or N (don't): ";

	 $response = <STDIN>;

	 chomp $response;
      }
      
      print "\n";
      
      if ( $response !~ /[nN]/ ) {
	 
	 # Grab the newest data from the server.
	 
	 &downloadNewRefSeqData();
	 
	 # Extract downloaded tarball into .genomeData.
	 
	 &expandNewRefSeqData_genomeData();

	 $newCopy = 1;

      } else {
	 
	 die("FATAL: Can't continue Phymm setup without genomic data; exiting.\n\n");
      }
   }
}

# Check for the existence of a gzipped tarball from the RefSeq FTP server,
# as well as any content in the .genomeData directory.

sub checkRefSeqData {
   
   if ( not -e '.genomeData' ) {
      
      # For whatever reason, the .genomeData directory doesn't (or has ceased to) exist.

      die("\nFATAL: Directory \".genomeData\" is missing; can't continue with setup.\n\n");

   } else {
      
      # The .genomeData directory exists, as expected.

      # If there's *anything at all* in the '.genomeData' directory, assume it's a genomic library and return 1.
   
      opendir DATA, ".genomeData";

      my @subs = readdir DATA;

      closedir DATA;

      foreach my $sub ( @subs ) {

	 if ( -d ".genomeData/$sub" and $sub !~ /^\./ ) {
	    
	    return 1;
	 }
      }
      
      # The .genomeData directory is empty.  Return 0.

      return 0;
   }
}

# Download a current tarball from RefSeq.

sub downloadNewRefSeqData {
   
   my $downloadLog = "$logDir/00_download_refseq_tarball__out.txt";
   
   open OUT, ">$downloadLog" or die("FATAL: Can't open $downloadLog for writing.\n");

   my $downloadErr = "$logDir/01_download_refseq_tarball__err.txt";

   open ERR, ">$downloadErr" or die("FATAL: Can't open $downloadErr for writing.\n");

   my $ftpSite = 'ftp.ncbi.nih.gov';

   print "Connecting to $ftpSite...";

   print OUT "Connecting to $ftpSite...";

   my $ftp = Net::FTP->new($ftpSite, Passive => 1, Debug => 0, Timeout => 2400) or die "FATAL: Can't connect to $ftpSite: $@";

   $ftp->login() or ( print ERR "Can't log in: " . $ftp->message . "\n" and close ERR and close OUT and die "FATAL: Can't log into $ftpSite: " . $ftp->message . "\n" );

   print "done.\n\n";

   print OUT "done.\n\n";

   my $targetDir = '/genomes/archive/old_refseq/Bacteria/';

   $ftp->cwd($targetDir) or ( print ERR "Can't change remote directory: " . $ftp->message . "\n" and close ERR and close OUT and die "FATAL: Can't change remote directory: " . $ftp->message . "\n" );

   my $targetFile = 'all.gbk.tar.gz';

   print "Downloading $targetFile (this will take a while; the file is over 2GB)...";

   print OUT "Downloading $targetFile (this will take a while; the file is over 2GB)...";

   $ftp->binary();
   
   $ftp->get($targetFile, "./$targetFile") or ( print ERR "Failed to download file: " . $ftp->message . "\n" and close ERR and close OUT and die "FATAL: Failed to download file: " . $ftp->message . "\n" );

   print "done.\n\n";

   print OUT "done.\n\n";

   $ftp->quit;

   close OUT;

   close ERR;
}

# Unzip/untar the downloaded genome-data file into the .genomeData directory.

sub expandNewRefSeqData_genomeData {
   
   my $expansionLog = "$logDir/02_NEW_refseq_tarball_expansion__out.txt";

   open OUT, ">$expansionLog" or die("FATAL: Can't open $expansionLog for writing.\n");

   my $expansionErr = "$logDir/03_NEW_refseq_tarball_expansion__err.txt";
   
   print "Expanding downloaded genome data into .genomeData/ for database overwrite (this'll take a while)...";

   print OUT "Expanding downloaded genome data into .genomeData/ for database overwrite:\n\n";

   print OUT "Running command: \"tar -zxvf all.gbk.tar.gz -C .genomeData > $expansionLog 2> $expansionErr\"\n\n";

   close OUT;

   system("tar -zxvf all.gbk.tar.gz -C .genomeData > $expansionLog 2> $expansionErr");

   open OUT, ">>$expansionLog" or die("FATAL: Can't open $expansionLog for appending.\n");
   
   print "done.\n\nRemoving compressed genome data...";
   
   print OUT "done.\n\nRemoving compressed genome data...";
   
   system("rm -rf all.gbk.tar.gz");
   
   print "done.\n\n";
   
   print OUT "done.\n\n";

   close OUT;
}

# Unzip/untar the downloaded genome-data file into a temp directory for comparison/update.

sub expandNewRefSeqData_tempDir {
   
   my $expansionLog = "$logDir/02_UPDATE_refseq_tarball_expansion__out.txt";

   open OUT, ">$expansionLog" or die("FATAL: Can't open $expansionLog for writing.\n");

   my $expansionErr = "$logDir/03_UPDATE_refseq_tarball_expansion__err.txt";
   
   print "Expanding downloaded genome data for database update (this'll take a while)...";
   
   print OUT "Expanding downloaded genome data for database update:\n\n";
   
   system("mkdir .tempGenomeData");

   print OUT "Running command: \"tar -zxvf all.gbk.tar.gz -C .tempGenomeData > $expansionLog 2> $expansionErr\"\n\n";

   close OUT;

   system("tar -zxvf all.gbk.tar.gz -C .tempGenomeData > $expansionLog 2> $expansionErr");

   open OUT, ">>$expansionLog" or die("FATAL: Can't open $expansionLog for appending.\n");
   
   print "done.\n\nRemoving tarball...";
   
   print OUT "\ndone.\n\nRemoving tarball...";
   
   system("rm -rf all.gbk.tar.gz");
   
   print "done.\n\n";
   
   print OUT "done.\n\n";

   close OUT;
}

# Standardize the names of the RefSeq organism directories to match the names given in the .gbk files.

sub renameNewRefSeqDirs {
   
   # Get the name of the genome data directory to scan.

   my $dataDir = shift;

   # Load known GenBank SOURCE-field inconsistencies for later processing.

   my $GBerrors = {};

   &loadGBerrors($GBerrors);

   my $warningsIssued = 0;

   my $renameLog = "$logDir/04_rename_genome_directories__out.txt";

   if ( $dataDir eq '.tempGenomeData' ) {
      
      $renameLog =~ s/04/06/;

      $renameLog =~ s/_genome/_temp_genome/;
   }

   open OUT, ">$renameLog" or die("FATAL: Can't open $renameLog for writing.\n");
   
   my $renameErr = "$logDir/05_rename_genome_directories__err.txt";

   if ( $dataDir eq '.tempGenomeData' ) {
      
      $renameErr =~ s/05/07/;

      $renameErr =~ s/_genome/_temp_genome/;
   }

   system("rm -f $renameErr");

   print "Standardizing organism directory names and data in $dataDir...";
   
   print OUT "Standardizing organism directory names and data in $dataDir...\n\n";
   
   opendir DATA, "$dataDir" or die("FATAL: Can't open $dataDir for scanning.\n");

   my @subs = grep { !/^\./ } sort { $a cmp $b } readdir DATA;

   closedir DATA;
   
   # Iterate over subdirectories of .genomeData.

   foreach my $sub ( @subs ) {
      
      # Organism-wide variable.  Used for sanity checking.

      my $orgName = '';
      
      if ( -d "$dataDir/$sub" ) {
	 
	 opendir SUB, "$dataDir/$sub" or die("FATAL: Can't open $dataDir/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;
	 
	 # Delete unwanted files.

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\%$/ ) {
	       
	       # This file's been obsolesced by NCBI, but for whatever reason remained in the distribution.  Kill it.
	       
	       system("rm $dataDir/$sub/$file");
	    }

	    # We only want accessions beginning with NC_ - all others are in varying stages of unreadiness.
	    
	    # Also, {NC_008265, NC_017594, NC_018285, NC_002179, NC_002180 and NC_017732} are phages and not bacteria; delete them.

	    if ( $file =~ /\.gbk$/ and ( $file !~ /NC_/ or $file =~ /NC_008265/ or $file =~ /NC_017594/
                                      or $file =~ /NC_018285/ or $file =~ /NC_002179/ or $file =~ /NC_002180/
                                      or $file =~ /NC_017732/ ) ) {
	       
	       system("rm $dataDir/$sub/$file");
	       
	    }

	 } # end ( loop deleting unwanted files )
	 
	 # Refresh the file list.  We may have deleted some.

	 opendir SUB, "$dataDir/$sub" or die("FATAL: Can't open $dataDir/$sub for scanning.\n");

	 @files = readdir SUB;

	 closedir SUB;
	 
	 # Go through each GenBank file in the current subdirectory of .genomeData.
	 # 
	 # 1. Grab the (possibly multi-line) organism name from the SOURCE field.
	 # 
	 # 2. If $orgName hasn't been set yet, set it to the current (retrieved) organism name.
	 # 
	 # 3. If $orgName has been set, then compare it to the current (retrieved) organism name.
	 #    If they don't match, it means there are different listings in the SOURCE field
	 #    across different files in the same subdirectory.  Check against a list of known
	 #    GenBank inconsistencies and fix if possible.  Otherwise, complain with a warning.
	 
	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\.gbk$/ ) {
	       
	       my $fullFile = "$dataDir/$sub/$file";
	       
	       open IN, "<$fullFile" or die("FATAL: Can't open $fullFile for reading.\n");
	       
	       my $found = 0;
	       
	       my $recording = 0;
	       
	       my $currentOrgName = '';
	       
	       while ( not $found and my $line = <IN> ) {
		  
		  chomp $line;

                  $line =~ s/\s*$//;

		  if ( $recording and $line !~ /^\s*ORGANISM/ ) {
		     
		     $line =~ /^\s*(.*)/;

		     $currentOrgName .= ' ' . $1;

		  } elsif ( $line =~ /^SOURCE\s*(.*)/ ) {
		     
		     $currentOrgName = $1;

		     $recording = 1;
		     
		  } elsif ( $line =~ /^\s*ORGANISM/ ) {
		     
		     $found = 1;
		  }
	       }

	       close IN;
	       
	       # Sanity check on organism name.
	       
	       $currentOrgName = &fixOrgName($currentOrgName, $GBerrors);

	       if ( $orgName eq '' ) {
		  
		  $orgName = $currentOrgName;

	       } elsif ( $orgName ne $currentOrgName ) {
                  
                  if ( $currentOrgName =~ /^$orgName/ ) {
                     
                     # Do nothing; we want to keep the version that's a prefix of the other.

                  } elsif ( $orgName =~ /^$currentOrgName/ ) {
                     
                     # Save the shorter name.

                     $orgName = $currentOrgName;

                  } else {
                     
                     # We really don't know how to resolve this.

                     open ERR, ">>$renameErr" or die("FATAL: Can't open $renameErr for appending.\n");

                     print ERR "WARNING: Inconsistent source-organism information listed across multiple .gbk files in genome directory \"$dataDir/$sub\": \"$orgName\" and \"$currentOrgName\".\n\n";

                     $warningsIssued = 1;

                     close ERR;
                  }
	       }
	       
	    } # end ( filename check for /\.gbk$/ )
	    
	 } # end ( file iterator on current subdirectory )
	 
	 # Check to ensure we haven't just emptied the current subdirectory.  If we have, delete it.

	 opendir SUB, "$dataDir/$sub" or die("FATAL: Can't open $dataDir/$sub for scanning.\n");

	 @files = readdir SUB;

	 closedir SUB;

	 my $found = 0;

	 foreach my $file ( @files ) {
	    
	    if ( $file ne '.' and $file ne '..' ) {
	       
	       $found = 1;
	    }
	 }

	 if ( not $found ) {
	    
	    # No files left in this subdirectory.  Delete it.
	    
	    system("rm -rf $dataDir/$sub\n");
	 }
	 
	 if ( $found ) {
	    
	    # The subdirectory isn't empty.  Continue.

	    # If we're not dealing with a completely alien directory - i.e., if we
	    # successfully found an $orgName - process the directory, otherwise ignore.

	    if ( $orgName ne '' ) {
	       
	       # Rename the current subdirectory to a standard name based on the species name quoted in its .gbk files.

	       my $newSubName = &toDirName($orgName);
	       
	       if ( $newSubName ne $sub ) {
		  
		  # We've changed the directory name.

		  if ( not -e "$dataDir/$newSubName" ) {
	       
		     # There's not yet a directory labeled with the new name.  Just move the old one over.
	       
		     my $tempSub = $sub;

                     $tempSub = &unixEscape($tempSub);

		     my $command = "mv $dataDir/$tempSub $dataDir/$newSubName\n";

		     system($command);

                     print OUT "   ...moved files from \"$sub\" to \"$newSubName\"...\n";

		  } else {
	       
		     # We've made a change, but the target directory (labeled with the
		     # new name) already exists; move the /files/ from the old dir to
		     # the new one.
	       
		     my $tempSub = $sub;

                     $tempSub = &unixEscape($tempSub);

		     my $command = "mv $dataDir/$tempSub/* $dataDir/$newSubName\n";

		     system($command);
	       
		     # Delete the now-empty directory labeled with the old name.
	       
		     $command = "rm -rf $dataDir/$tempSub\n";

		     system($command);

                     print OUT "   ...moved files from \"$sub\" to \"$newSubName\"...\n";
		  }

	       } # end ( check to see if we changed the directory name from what it was )

	    } # end ( check to discover whether or not we've been able to establish an $orgName )
	    
	 } # end ( check to ensure current subdirectory isn't empty )
	 
      } # end ( filter on /^\./ )
      
   } # end foreach ( subdirectory in .genomeData )
   
   print OUT "\ndone.\n\n";

   close OUT;
   
   print "done.";

   if ( $warningsIssued ) {
      
      print "\n\n   NOTE: non-fatal warnings were issued while standardizing names of\n";
      
      print "   individual organisms' genome data directories; this is almost always\n";
      
      print "   caused by a typo or syntax inconsistency in GenBank-derived metadata.\n";
      
      print "   Check the file \"$renameErr\" after\n";
      
      print "   setup completes to see details of these warnings.";
   }

   print "\n\n";
}

# Generate all needed metadata for later Phymm processing.

sub generateRefSeqMetafiles {
   
   my $taxonomyLog = "$logDir/06_generate_taxonomic_metadata__out.txt";

   if ( $newCopy == 2 ) {
      
      $taxonomyLog =~ s/06/10/;
   }

   open OUT, ">$taxonomyLog" or die("FATAL: Can't open $taxonomyLog for writing.\n");

   my $taxonomyErr = "$logDir/07_generate_taxonomic_metadata__err.txt";

   if ( $newCopy == 2 ) {
      
      $taxonomyErr =~ s/07/11/;
   }

   system("rm -f $taxonomyErr");

   print "Generating taxonomic metafiles...";
   
   print OUT "Generating taxonomic metafiles...\n\n";
   
   my $species = {};
   
   my $definition = {};
   
   my $dataDir = '.genomeData';

   opendir DATA, "$dataDir" or die("FATAL: Can't open $dataDir for scanning.\n");

   my @subs = grep { !/^\./ } sort { $a cmp $b } readdir DATA;

   closedir DATA;
   
   # Go through each subdirectory of .genomeData and grab the defline.  Reconstruct the organism name from the directory name.
   
   print OUT "   Grabbing deflines from .gbk files...";

   foreach my $sub ( @subs ) {
      
      my $defLine = '';
      
      my $orgName = '';
      
      if ( -d "$dataDir/$sub" ) {
	 
	 opendir SUB, "$dataDir/$sub" or die("FATAL: Can't open $dataDir/$sub for scanning.\n");

	 my @files = readdir SUB;

	 closedir SUB;
	 
	 $orgName = &fromDirName($sub);

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /(.+)\.gbk$/ ) {
	       
	       my $prefix = $1;
	       
	       my $fullFile = "$dataDir/$sub/$file";

	       open IN, "<$fullFile" or die("FATAL: Can't open $fullFile for reading.\n");
	       
	       my $found = 0;
	       
	       $defLine = '';
	       
	       my $recording = 0;
	       
	       while ( not $found and my $line = <IN> ) {
		  
		  chomp $line;

                  $line =~ s/\s*$//;

		  if ( $recording and $line !~ /^ACCESSION/ ) {
		     
		     $line =~ /^\s*(.*)/;

		     $defLine .= ' ' . $1;

		  } elsif ( $line =~ /^DEFINITION\s*(.*)/ ) {
		     
		     $defLine = $1;

		     $recording = 1;
		     
		  } elsif ( $line =~ /^ACCESSION/ ) {
		     
		     $found = 1;
		  }
	       }
	       
	       $species->{$prefix} = $orgName;

	       $definition->{$prefix} = $defLine;
	       
	    } # end ( filename check for /\.gbk$/ )
	    
	 } # end ( file iterator on current subdirectory )
	 
      } # end ( filter on /^\./ )
      
   } # end foreach ( subdirectory in .genomeData )

   print OUT "done.\n";
   
   # Output all the saved data to an accession map file.

   my $idMapFile = '.taxonomyData/.0_accessionMap/accessionMap.txt';

   open IDMAP, ">$idMapFile" or die("FATAL: Can't open $idMapFile for writing.\n");

   print OUT "   Writing accession-to-species map to \"$idMapFile\"...";
   
   foreach my $prefix ( sort { $species->{$a} cmp $species->{$b} || $definition->{$a} cmp $definition->{$b} } keys %$species ) {
      
      my $orgName = $species->{$prefix};

      my $defLine = $definition->{$prefix};

      my $seqType = '';

      if ( $defLine =~ /chromosom/ or $defLine =~ /complete\s+genome/ ) {
	 
	 $seqType = 'chromosome';

      } elsif ( $defLine =~ /plasmid/ ) {
	 
	 $seqType = 'plasmid';

      } elsif ( $defLine =~ /complete\s+sequence/ ) {
	 
	 $seqType = 'chromosome';

      } else {
	 
	 $seqType = 'chromosome';
      }

      my $newSubName = &toDirName($orgName);

      print IDMAP "$newSubName\t$prefix\t$seqType\t$defLine\n";
   }
   
   close IDMAP;

   print OUT "done.\n";
   
   # Grab the NCBI taxonomy, parse, and save in a format we can cope with.

   print OUT "   Scraping NCBI taxonomy webpages for bacteria and archaea into .taxonomyData/.1_fromNCBI/ ...";

   system(".scripts/.taxonomyInfo/0_getNCBIpages.pl");

   print OUT "done.\n   Parsing HTML (bacteria)...";
   
   system(".scripts/.taxonomyInfo/1_stripHTML_bacteria.pl 2>> $taxonomyErr");

   print OUT "done.\n   Extracting taxonomic tree (bacteria)...";
   
   system(".scripts/.taxonomyInfo/2_getTree_bacteria.pl 2>> $taxonomyErr");

   print OUT "done.\n   Parsing HTML (archaea)...";
   
   system(".scripts/.taxonomyInfo/3_stripHTML_archaea.pl 2>> $taxonomyErr");

   print OUT "done.\n   Extracting taxonomic tree (archaea)...";

   system(".scripts/.taxonomyInfo/4_getTree_archaea.pl 2>> $taxonomyErr");

   print OUT "done.\n   Creating flat taxonomy files...";

   system(".scripts/.taxonomyInfo/5_configureTaxonomicData.pl 2>> $taxonomyErr");

   print OUT "done.\n";
   
   print "done.\n\n";

   print OUT "\n...taxonomic metadata generation complete.\n";

   close OUT;
}

# Convert the .gbk files that came with the RefSeq distribution to corresponding .fna files which contain only the nucleotide sequence from the .gbk file.

sub genbankToFasta {
   
   my $fastaLog = "$logDir/08_genbank_to_fasta__out.txt";

   if ( $newCopy == 2 ) {
      
      $fastaLog =~ s/08/12/;
   }

   my $fastaErr = "$logDir/09_genbank_to_fasta__err.txt";

   if ( $newCopy == 2 ) {
      
      $fastaErr =~ s/09/13/;
   }

   open OUT, ">$fastaLog" or die("FATAL: Can't open $fastaLog for writing.\n");

   print "Extracting nucleotide sequences from .gbk files into .fna files as necessary...";
      
   print OUT "Extracting nucleotide sequences from .gbk files into .fna files as necessary...\n\n";

   close OUT;

   system(".scripts/genbankToFasta.pl .genomeData >> $fastaLog 2> $fastaErr");

   print "done.\n\n";

   open OUT, ">>$fastaLog" or die("FATAL: Can't open $fastaLog for appending.\n");

   print OUT "\ndone.\n";

   close OUT;
}

# Create training files for the ICMs if they don't already exist.

sub createTrainingFiles {
   
   my $trainFileLog = "$logDir/10_create_training_files__out.txt";

   if ( $newCopy == 2 ) {
      
      $trainFileLog =~ s/10/14/;
   }

   open OUT, ">$trainFileLog" or die("FATAL: Can't open $trainFileLog for writing.\n");

   my $trainFileErr = "$logDir/11_create_training_files__err.txt";

   if ( $newCopy == 2 ) {
      
      $trainFileErr =~ s/11/15/;
   }

   print "Creating training files for ICMs as necessary...";
   
   print OUT "Creating training files for ICMs as necessary...\n\n";

   close OUT;

   system(".scripts/createTrainingFiles.pl .genomeData 500 >> $trainFileLog 2> $trainFileErr");

   print "done.\n\n";

   open OUT, ">>$trainFileLog" or die("FATAL: Can't open $trainFileLog for appending.\n");

   print OUT "\ndone.\n";

   close OUT;
}

# Build the ICMs which form the core of Phymm's classification system

sub buildICMs {
   
   my $buildICMsLog = "$logDir/14_build_ICMs__out.txt";

   if ( $newCopy == 2 ) {
      
      $buildICMsLog =~ s/14/18/;
   }

   open OUT, ">$buildICMsLog" or die("FATAL: Can't open $buildICMsLog for writing.\n");

   my $buildICMsErr = "$logDir/15_build_ICMs__err.txt";

   if ( $newCopy == 2 ) {
      
      $buildICMsErr =~ s/15/19/;
   }

   system("rm -f $buildICMsErr");

   # First, check to see if the ICM code from the Glimmer package has been untarred and/or compiled yet.  Do whatever needs doing.

   if ( not -e '.scripts/.icmCode/LICENSE' ) {
      
      # The package hasn't been uncompressed yet.  Do so.
      
      print "Uncompressing ICM code...";
      
      print OUT "Uncompressing ICM code...\n\n";

      close OUT;
      
      system("tar -zxvf .scripts/.icmCode/glimmer3_plusSS.tgz -C .scripts/.icmCode >> $buildICMsLog 2>> $buildICMsErr");

      open OUT, ">>$buildICMsLog" or die("FATAL: Can't open $buildICMsLog for appending.\n");
      
      print "done.\n\n";

      print OUT "\ndone.\n\n";
   }
   
   # Next, resolve a conflict based on mutually exclusive syntax requirements across the gcc 4.3/4.4 barrier.

   open (CMD, "(gcc -v) 2>&1 |");

   while ( my $line = <CMD> ) {
   
      if ( $line =~ /gcc\s+version\s+(\S+)/ ) {
      
	 my $versionString = $1;
      
	 my @versionNums = ();

	 while ( $versionString =~ /^(\d+)\./ ) {
	 
	    my $nextNum = $1;

	    push @versionNums, $nextNum;

	    $versionString =~ s/^$nextNum\.//;
	 }
      
	 if ( $versionString =~ /^(\d+)/ ) {
	 
	    push @versionNums, $1;
	 }

	 my $topNum = $versionNums[0];
      
	 my $secondNum = -1;
      
	 if ( $#versionNums > 0 ) {
	 
	    $secondNum = $versionNums[1];
	 }

	 my $newGCC = 0;

	 if ( $topNum > 4 or ( $topNum == 4 and $secondNum >= 4 ) ) {
	 
	    $newGCC = 1;
	 }
	 
	 # Rename the appropriate files, if we haven't already done this.

	 if ( $newGCC == 0 ) {
            
            print OUT "Using source files for GCC v[<=4.3.*].\n";
	    
	    # We've got an older version of gcc.

	    if ( not -e '.scripts/.icmCode/SimpleMake/icm.cc' ) {
	       
	       system("mv .scripts/.icmCode/SimpleMake/icm_oldgcc.cc .scripts/.icmCode/SimpleMake/icm.cc");
	       system("mv .scripts/.icmCode/SimpleMake/gene_oldgcc.cc .scripts/.icmCode/SimpleMake/gene.cc");
	       system("mv .scripts/.icmCode/SimpleMake/delcher_oldgcc.hh .scripts/.icmCode/SimpleMake/delcher.hh");
	    }

	 } else {
	    
            print OUT "Using source files for GCC v[>=4.4.*].\n";
	    
	    # We've got a newer version of gcc.

	    if ( not -e '.scripts/.icmCode/SimpleMake/icm.cc' ) {
	       
	       system("mv .scripts/.icmCode/SimpleMake/icm_newgcc.cc .scripts/.icmCode/SimpleMake/icm.cc");
	       system("mv .scripts/.icmCode/SimpleMake/gene_newgcc.cc .scripts/.icmCode/SimpleMake/gene.cc");
	       system("mv .scripts/.icmCode/SimpleMake/delcher_newgcc.hh .scripts/.icmCode/SimpleMake/delcher.hh");
	    }
	 }
      }
   }

   close CMD;

   # Now we can compile.

   if ( not -e ".scripts/.icmCode/bin/build-icm" ) {
      
      # The code hasn't been compiled yet.  Do so.
      
      print "Compiling ICM code...";
      
      print OUT "Compiling ICM code...\n\n";

      close OUT;
      
      chdir(".scripts/.icmCode/SimpleMake");

      system("make >> ../../../$buildICMsLog 2>> ../../../$buildICMsErr");

      chdir("../../..");

      open OUT, ">>$buildICMsLog" or die("FATAL: Can't open $buildICMsLog for appending.\n");
      
      print "done.\n\n";
      
      print OUT "\ndone.\n\n";
   }
   
   # Sanity: check for one of the more critical binaries.

   if ( not -e ".scripts/.icmCode/bin/build-icm" ) {
      
      die("FATAL: The ICM code does not appear to have compiled successfully.\n\nPlease try to compile manually (cd .scripts/.icmCode/SimpleMake; make);\nif that doesn't work, check to make sure you have\nthe appropriate version of gcc (gcc -v) for the version of Phymm\nthat you've downloaded (see website for details).\n");
   }

   # If ICMs haven't yet been built, build them.

   close OUT;

   system(".scripts/buildICMs.pl .genomeData $buildICMsLog $buildICMsErr");
}

# Build a local BLAST database from the RefSeq FASTA files.

sub generateBlastDB {
   
   my $locateBlastBinary = `which makeblastdb 2>&1`;
   
   if ( $locateBlastBinary =~ /no\s+makeblastdb\s+in/ or $locateBlastBinary =~ /command\s+not\s+found/i ) {
      
      die("FATAL: You must have a local copy of the BLAST software installed and accessible via your \"PATH\" environment variable.  Please see the README file for details.\n\n");

   } else {
      
      # Standalone BLAST is installed.  Do we need to regenerate the BLAST DB?  If there's no DB, generate one.
      
      my $regenDB = 0;

      if ( not -e ".blastData/phymm_BLAST_DB.nal" ) {
	 
	 $regenDB = 1;

      } else {
	 
	 # If there is one, search the lastmod timestamps of the .fna files in .genomeData.
	 # If we find one that's newer than the blast DB, regenerate the DB.
	 # 
	 # We won't check .fna files in .genomeData/.userAdded here, because if the user has
	 # used addCustomGenome.pl, then they've been repeatedly instructed to remember to
	 # rebuild the BLAST DB when they were done (if they selected the option to do it by
	 # hand in the first place, that is; by default, it's rebuilt automatically).  If they
	 # decided not to do it despite repeated cautions, well then who knows, they might
	 # have their reasons.

	 my $blastTimeStamp = (stat('.blastData/phymm_BLAST_DB.nal'))[9];
	 
	 my $genomeDir = '.genomeData';

	 opendir DOT, $genomeDir or die("FATAL: Can't open $genomeDir for scanning.\n");

	 my @subs = grep { !/^\./ } readdir DOT;

	 closedir DOT;

	 OUTER: foreach my $sub ( @subs ) {
	    
	    my $fullSub = "$genomeDir/$sub";

	    if ( -d $fullSub ) {
	       
	       opendir DOT, $fullSub or die("FATAL: Can't open $fullSub for scanning.\n");

	       my @files = grep { /\.fna$/ } readdir DOT;

	       closedir DOT;

	       foreach my $file ( @files ) {
		  
		  my $fullFile = "$fullSub/$file";

		  my $currentTimeStamp = (stat($fullFile))[9];

		  if ( $currentTimeStamp > $blastTimeStamp ) {
			
		     # Found one.  Stop looking.

		     $regenDB = 1;

		     last OUTER;
		  }
	       }
	    }
	 }
      }

      if ( $regenDB ) {
	 
         my $blastLog = "$logDir/12_generate_BLAST_DB__out.txt";

         if ( $newCopy == 2 ) {
            
            $blastLog =~ s/12/16/;
         }

         open LOG, ">$blastLog" or die("FATAL: Can't open $blastLog for writing.\n");

         my $blastErr = "$logDir/13_generate_BLAST_DB__err.txt";

         if ( $newCopy == 2 ) {
            
            $blastErr =~ s/13/17/;
         }

         system("rm -f $blastErr");

	 # Concatenate all genome files and build a local BLAST DB from them.
	 
	 print "Generating local BLAST database from genome data...";
      
	 print LOG "Generating local BLAST database from genome data...\n\n";
      
	 # First, grab a list of relative paths to all the RefSeq genomic FASTA files.
      
	 my @fastaFiles = ();
      
	 opendir DOT, '.genomeData' or die("FATAL: Can't open .genomeData for scanning.\n");

	 my @subs = grep { !/^\./ } sort { $a cmp $b } readdir DOT;

	 closedir DOT;

	 foreach my $sub ( @subs ) {
	 
	    if ( -d ".genomeData/$sub" ) {
	    
	       opendir DOT, ".genomeData/$sub" or die("FATAL: Can't open .genomeData/$sub for scanning.\n");

	       my @files = sort { $a cmp $b } grep { /\.fna$/ } readdir DOT;

	       closedir DOT;

	       foreach my $file ( @files ) {
	       
		  push @fastaFiles, ".genomeData/$sub/$file";
	       }
	    }
	 }
	 
	 # Next, grab a list of relative paths to any user-added FASTA files.
	 
	 my $userDir = '.genomeData/.userAdded';
	 
	 if ( -e $userDir ) {
	    
	    opendir DOT, $userDir or die("FATAL: Can't open $userDir for scanning.\n");

	    @subs = grep { !/^\./ } sort { $a cmp $b } readdir DOT;

	    closedir DOT;

	    foreach my $sub ( @subs ) {
	       
	       if ( -d "$userDir/$sub" ) {
		  
		  opendir DOT, "$userDir/$sub" or die("FATAL: Can't open $userDir/$sub for scanning.\n");

		  my @files = sort { $a cmp $b } grep { /\.fna$/ } readdir DOT;

		  closedir DOT;

		  foreach my $file ( @files ) {
		     
		     push @fastaFiles, "$userDir/$sub/$file";
		  }
	       }
	    }
	 }

	 # Concatenate all the FASTA files together for database creation.
      
	 my $command = 'cat ';

	 my $argCount = 0;
	 
	 my $wrapped = 0;

         my $fileIndex = -1;

	 foreach my $file ( @fastaFiles ) {
	    
            $fileIndex++;

	    $command .= "$file \\\n";

	    $argCount++;

	    if ( $argCount % 10 == 0 and $fileIndex != $#fastaFiles ) {
                  
               if ( not $wrapped ) {
                     
                  $wrapped = 1;

                  $command .= "> blastDB_data.txt\n";

               } else {
                     
                  $command .= ">> blastDB_data.txt\n";
               }
                  
               $command .= 'cat ';
	    }
	 }
	 
	 if ( not $wrapped ) {
	    
	    $command .= "> blastDB_data.txt\n";

	 } else {
	    
	    $command .= ">> blastDB_data.txt\n";
	 }
	 
	 # There's a problem with using system() for a command this long, apparently, so we use a hack.

	 my $shLoc = `which sh`;

	 chomp $shLoc;

	 my $outFile = "command_temp.sh";

	 open OUT, ">$outFile" or die("FATAL: Can't open $outFile for writing.\n");

	 print OUT "#!$shLoc\n\n$command";

	 close OUT;

	 system("chmod 775 command_temp.sh");

	 system("./command_temp.sh");

	 system("rm command_temp.sh");

	 # Create the local database copy.

         close LOG;
	 
	 system("makeblastdb -in blastDB_data.txt -title phymm_BLAST_DB -out phymm_BLAST_DB -dbtype nucl >> $blastLog 2> $blastErr");

	 system("mv phymm_BLAST_DB.* .blastData");

	 system("rm blastDB_data.txt");
      
	 print "done.\n\n";

         open LOG, ">>$blastLog" or die("FATAL: Can't open $blastLog for appending.\n");

         print LOG "\ndone.\n\n";

         close LOG;

      } # end if ( we're generating a DB at all )

   } # end if ( local copy of BLAST binary is installed )
}

sub toDirName {
   
   my $name = shift;

   my $result = $name;

   $result =~ s/\_/UNDERSCORE/g;

   $result =~ s/\+/PLUS/g;

   $result =~ s/\s/\_/g;

   $result =~ s/\//SLASH/g;

   $result =~ s/\(/LPAREN/g;

   $result =~ s/\)/RPAREN/g;

   $result =~ s/\'/SINGLEQUOTE/g;

   $result =~ s/\"/DOUBLEQUOTE/g;
	    
   $result =~ s/\:/COLONCOLON/g;

   $result =~ s/\;/SEMICOLON/g;

   return $result;
}

sub fromDirName {
   
   my $name = shift;

   my $result = $name;

   $result =~ s/\_/ /g;

   $result =~ s/UNDERSCORE/\_/g;

   $result =~ s/PLUS/\+/g;

   $result =~ s/SLASH/\//g;

   $result =~ s/LPAREN/\(/g;

   $result =~ s/RPAREN/\)/g;

   $result =~ s/SINGLEQUOTE/\'/g;

   $result =~ s/DOUBLEQUOTE/\"/g;
	    
   $result =~ s/COLONCOLON/\:/g;

   $result =~ s/SEMICOLON/\;/g;

   return $result;
}

sub loadGBerrors {
   
   my $errors = shift;

   my $errorFile = '.taxonomyData/knownGBerrors.txt';

   open IN, "<$errorFile" or die("FATAL: Can't open $errorFile for reading.\n");

   while ( my $line = <IN> ) {
      
      chomp $line;

      if ( $line =~ /\t/ ) {
         
         (my $canonName, my $pattern) = split(/\t/, $line);

         $errors->{$pattern} = $canonName;
      }
   }

   close IN;
}

sub fixOrgName {
   
   my $orgName = shift;

   my $errors = shift;

   foreach my $pattern ( keys %$errors ) {
      
      my $correctName = $errors->{$pattern};

      if ( $orgName =~ /$pattern/ ) {
	 
	 $orgName = $correctName;

	 return $orgName;
      }
   }

   return $orgName;
}

sub updateGenomeData {
   
   # Classify each genome directory in the new RefSeq set as one of the following:
   #    [new organism added since existing installation was last updated],
   #    [organism contained in existing installation, but /hasn't/ been changed], or
   #    [organism contained in existing installation and /has/ been changed].

   my $updateLog = "$logDir/08_update_genome_data__out.txt";

   open OUT, ">$updateLog" or die("FATAL: Can't open $updateLog for writing.\n");

   my $updateErr = "$logDir/09_update_genome_data__err.txt";

   open ERR, ">$updateErr" or die("FATAL: Can't open $updateErr for writing.\n");

   # Load the existing genome directory list, along with .fna file sizes.

   my $extantGenomes = {};

   opendir GENOMEDATA, '.genomeData' or die("FATAL: Can't open .genomeData for scanning.\n");

   my @subs = readdir GENOMEDATA;

   closedir GENOMEDATA;

   foreach my $sub ( @subs ) {
      
      if ( -d ".genomeData/$sub" and $sub !~ /^\./ ) {
	 
	 opendir DOT,  ".genomeData/$sub" or die("FATAL: Can't open .genomeData/$sub for scanning.\n");

	 my @files = readdir DOT;

	 closedir DOT;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\.fna$/ ) {
	       
	       my $fileSize = -s ".genomeData/$sub/$file";

	       $extantGenomes->{$sub}->{$file} = $fileSize;
	    }
	 }
      }
   }
   
   # Create .fna files for the new download.

   print "Creating temporary FASTA files from new .gbk files for update check...";

   print OUT "Creating temporary FASTA files from new .gbk files for update check...\n\n";

   close OUT;

   close ERR;

   system(".scripts/genbankToFasta.pl .tempGenomeData >> $updateLog 2>> $updateErr");

   open OUT, ">>$updateLog" or die("FATAL: Can't open $updateLog for appending.\n");

   open ERR, ">>$updateErr" or die("FATAL: Can't open $updateErr for appending.\n");

   print "done.\n\n";

   print OUT "done.\n\n";

   # Load the new genome directory list, along with .fna file sizes.

   my $newGenomes = {};

   opendir GENOMEDATA, '.tempGenomeData' or die("FATAL: Can't open .tempGenomeData for scanning.\n");

   my @subs = readdir GENOMEDATA;

   closedir GENOMEDATA;

   foreach my $sub ( @subs ) {
      
      if ( -d ".tempGenomeData/$sub" and $sub !~ /^\./ ) {
	 
	 opendir DOT,  ".tempGenomeData/$sub" or die("FATAL: Can't open .tempGenomeData/$sub for scanning.\n");

	 my @files = readdir DOT;

	 closedir DOT;

	 foreach my $file ( @files ) {
	    
	    if ( $file =~ /\.fna$/ ) {
	       
	       my $fileSize = -s ".tempGenomeData/$sub/$file";

	       $newGenomes->{$sub}->{$file} = $fileSize;
	    }
	 }
      }
   }

   # First, find the completely new critters, if any.
   
   my $addedGenomeHash = {};

   foreach my $newGenome ( keys %$newGenomes ) {
      
      if ( not exists($extantGenomes->{$newGenome}) ) {
	 
	 $addedGenomeHash->{$newGenome} = 1;
      }
   }
   
   my @addedGenomes = sort { $a cmp $b } keys %$addedGenomeHash;

   if ( $#addedGenomes == -1 ) {
      
      print "No new genomes to add.\n\n";

      print OUT "No new genomes to add.\n\n";

   } else {
      
      print "Adding new genomes:\n\n";

      print OUT "Adding new genomes:\n\n";

      foreach my $newGenome ( @addedGenomes ) {
	 
	 print "   $newGenome\n";
	 
	 print OUT "   $newGenome\n";
      }

      print "\n";

      print OUT "\n";
   }

   # Next, identify existing genomes which need to be updated, if any.
   
   my $updatedGenomeHash = {};

   foreach my $newGenome ( keys %$newGenomes ) {
      
      foreach my $file ( keys %{$newGenomes->{$newGenome}} ) {
	 
	 if ( not exists($extantGenomes->{$newGenome})
	   or not exists($extantGenomes->{$newGenome}->{$file})
	   or $extantGenomes->{$newGenome}->{$file} != $newGenomes->{$newGenome}->{$file} ) {
	    
	    if ( not exists($addedGenomeHash->{$newGenome}) ) {
	       
	       $updatedGenomeHash->{$newGenome}->{$file} = 1;
	    }
	 }
      }
   }

   my @updatedGenomes = sort { $a cmp $b } keys %$updatedGenomeHash;

   if ( $#updatedGenomes == -1 ) {
      
      print "No updates found for existing genomes.\n\n";

      print OUT "No updates found for existing genomes.\n\n";

   } else {
      
      print "Updating existing genomes with new data:\n\n";

      print OUT "Updating existing genomes with new data:\n\n";

      foreach my $updateGenome ( @updatedGenomes ) {
	 
	 print "   $updateGenome\n";
	 
	 print OUT "   $updateGenome\n";
      }

      print "\n";

      print OUT "\n";
   }

   # Process the additions: just move the directories over.

   print OUT "Processing new genomes...";

   foreach my $addGenome ( @addedGenomes ) {
      
      system("cp -R .tempGenomeData/$addGenome .genomeData");
   }

   print OUT "done.\n";
   
   # Process the updates: remove old .gbk files, as well
   # as any files generated from them, and copy the new
   # .gbk files over.

   print OUT "Processing existing-genome updates...";
   
   foreach my $updateGenome ( @updatedGenomes ) {
      
      my @updateFiles = keys %{$updatedGenomeHash->{$updateGenome}};

      foreach my $updateFile ( @updateFiles ) {
	 
	 $updateFile =~ /^(.+)\.fna/;

	 my $prefix = $1;

	 my $command = "rm -rf .genomeData/$updateGenome/$prefix\*";

	 system($command);

	 $command = "cp .tempGenomeData/$updateGenome/$prefix\.gbk .genomeData/$updateGenome/$prefix\.gbk";

	 system($command);
      }
   }

   print OUT "done.\n";

   # Delete the temporary genome repository; we're done with it.
   
   print OUT "Deleting temporary genome repository...";

   system("rm -rf .tempGenomeData");

   print OUT "done.\n";

   close OUT;

   close ERR;
}

sub unixEscape {
   
   my $arg = shift;

   $arg =~ s/\(/\\\(/g;

   $arg =~ s/\)/\\\)/g;

   $arg =~ s/\'/\\\'/g;

   $arg =~ s/\"/\\\"/g;

   $arg =~ s/ /\\ /g;

   return $arg;
}
