use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Statistics::FisherPitman',
    VERSION_FROM      => 'lib/Statistics/FisherPitman.pm', # finds $VERSION
    PREREQ_PM         => {
	List::Util => 1.18,
	Statistics::Descriptive => 2.6,
	},
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT_FROM  => 'lib/Statistics/FisherPitman.pm', # retrieve abstract from module
       AUTHOR         => 'R. Garton <rgarton AT cpan DOT org>') : ()),
);
