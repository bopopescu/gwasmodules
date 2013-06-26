#!/usr/bin/env python
"""
Examples:
	#2012.6.27, no clustering (-C 1), all downloaded files will be in gxfer3_74079526293008 subfolder (-F ..)
	# cut off 2 layers of folders from the URL (-c2) as the URL has two folders after its top-level URL.
	%s -I https://xfer.genome.wustl.edu/gxfer3/74079526293008/ -u aufiewisiuch -p uazeiraiquae
		-o /tmp/workflow.xml -F gxfer3_74079526293008 -j hcondor -l hcondor -C 1 -c2
	%s

2012.6.27
	a common class for pegasus workflows that work on NGS (next-gen sequencing) data

"""
import sys, os, math
__doc__ = __doc__%(sys.argv[0], sys.argv[0])

sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))

from pymodule import ProcessOptions, getListOutOfStr, PassingData, yh_pegasus, utils
from Pegasus.DAX3 import *
from AbstractWorkflow import AbstractWorkflow

class WgetWorkflow(AbstractWorkflow):
	__doc__ = __doc__
	option_default_dict = AbstractWorkflow.option_default_dict.copy()
	option_default_dict.update({
						("inputFname", 0, ): ["", 'i', 1, 'the input file that contains the URLs to wget, not supported yet.'],\
						("inputURL", 0, ): ["", 'I', 1, 'the URL and all the links (files only) embedded in the URL will be downloaded.'],\
						("pathToWget", 1, ): ['/usr/bin/wget', '', 1, 'path to the wget binary'],\
						("username", 0, ): ['', 'u', 1, 'username to access the URL'],\
						("password", 0, ): ['', 'p', 1, 'password to access the URL'],\
						("cut_dir_number", 1, int): [1, 'c', 1, 'how many levels of folders after the URL to be ignored. \
	i.e. Files in ftp://ftp.xemacs.org/pub/xemacs/ (cut_dir_number=1) will have prefix xemacs/. \
	It (--cut-dirs=..) is used in conjunction with -nH.'],\
						
						})
						#('bamListFname', 1, ): ['/tmp/bamFileList.txt', 'L', 1, 'The file contains path to each bam file, one file per line.'],\

	def __init__(self,  **keywords):
		"""
		2012.6.27
		"""
		AbstractWorkflow.__init__(self, **keywords)
	
	def addWgetJob(self, executable=None, url=None, relativePath=None, username=None, password=None,\
				targetFolder=None, logFile=None, cut_dir_number=1, parentJobLs=[], extraDependentInputLs=[], transferOutput=False, \
				extraArguments=None, job_max_memory=2000, **keywords):
		"""
		2012.6.27
		"""
		extraArgumentList = ['--user=%s'%(username), '--password=%s'%(password), '--recursive', '--no-parent',\
					'--continue', "--reject='index.html*'", "-nc -nH --cut-dirs=%s"%(cut_dir_number), "-P %s"%(targetFolder), \
					"%s/%s"%(url, relativePath)]
		
		"""
		# unlike -nd, --cut-dirs does not lose with subdirectories---for instance, with
		# -nH --cut-dirs=1, a beta/ subdirectory will be placed to xemacs/beta, as one would expect.
		
		-c
		--continue
		   Continue getting a partially-downloaded file.  This is useful when you want to finish up a download started
		   by a previous instance of Wget, or by another program. 
		
		-nc
		--no-clobber
		   If a file is downloaded more than once in the same directory, Wget's behavior depends on a few options,
		   including -nc.  In certain cases, the local file will be clobbered, or overwritten, upon repeated download.
		   In other cases it will be preserved.		
		
		-nd
		--no-directories
		   Do not create a hierarchy of directories when retrieving recursively.  With this option turned on, all
		   files will get saved to the current directory, without clobbering (if a name shows up more than once, the
		   filenames will get extensions .n).
		
		-np
		--no-parent
		   Do not ever ascend to the parent directory when retrieving recursively.  This is a useful option, since it
		   guarantees that only the files below a certain hierarchy will be downloaded.
	
		-nH
		--no-host-directories
		   Disable generation of host-prefixed directories.  By default, invoking Wget with -r http://fly.srk.fer.hr/
		   will create a structure of directories beginning with fly.srk.fer.hr/.  This option disables such behavior.

		-P prefix
		--directory-prefix=prefix
		   Set directory prefix to prefix.  The directory prefix is the directory where all other files and
		   subdirectories will be saved to, i.e. the top of the retrieval tree.  The default is . (the current
		   directory)
"""
		if extraArguments:
			extraArgumentList.append(extraArguments)
		#wget will add some portion of the URL path to the final output files depending on the cut_dir_number
		import urlparse
		url_path_list = urlparse.urlparse(url).path.split('/')[1:]	#[0] is empty because the path starts with '/' 
		subPath = '/'.join(url_path_list[cut_dir_number:])
		
		if relativePath.find('/')>=0:	#2013.06.26 it's a folder itself. so no straight output.
			sys.stderr.write("\n\t Warning: item %s is a folder, will not be staged out. You have to manually copy them out of scratch folder.\n"%\
							(relativePath))
			extraOutputLs = None
		else:
			extraOutputLs = [File(os.path.join(targetFolder, os.path.join(subPath, relativePath)))]
		#2012.6.27 don't pass the downloaded outputFile to argument outputFile of addGenericJob()
		# because it will add "-o" in front of it. "-o" of wget is reserved for logFile.
		return self.addGenericJob(executable=executable, inputFile=None, outputFile=logFile, \
						parentJobLs=parentJobLs, extraDependentInputLs=extraDependentInputLs, \
						extraOutputLs=extraOutputLs,\
						transferOutput=transferOutput, \
						extraArgumentList=extraArgumentList, job_max_memory=job_max_memory)
	
	def getFilenamesToBeDownloaded(self, url=None, username=None, password=None):
		"""
		2012.6.27
			get all links from within that url
		"""
		sys.stderr.write('Getting filenames from the url page ... ')
		
		# create a password manager
		import urllib2
		password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
		
		# Add the username and password.
		# If we knew the realm, we could use it instead of None.
		top_level_url = url
		password_mgr.add_password(None, top_level_url, username, password)
		
		handler = urllib2.HTTPBasicAuthHandler(password_mgr)
		
		# create "opener" (OpenerDirector instance)
		opener = urllib2.build_opener(handler)
		
		# use the opener to fetch a URL
		f = opener.open(url)
		
		"""
		# Install the opener.
		# Now all calls to urllib2.urlopen use our opener.
		urllib2.install_opener(opener)
		"""
		url_content = f.read()
		
		from HTMLParser import HTMLParser
		# create a subclass and override the handler methods
		class MyHTMLParser(HTMLParser):
							
			startTagList = []
			endTagList = []
			hrefValueList = []
			relativePathList = []
			def __init(self):
				"""
				2012.6.27 example http page source. lots of hrefs that are irrelevant.
					One rule is that the href attribute value should be equal to the data inside the <a> tag.
				
				<h1>Index of /gxfer3/74079526293008</h1>
<table><tr><th><img src="/icons/blank.gif" alt="[ICO]"></th><th><a href="?C=N;O=D">Name</a></th>
	<th><a href="?C=M;O=A">Last modified</a></th><th><a href="?C=S;O=A">Size</a></th>
	<th><a href="?C=D;O=A">Description</a></th></tr><tr><th colspan="5"><hr></th></tr>
<tr><td valign="top"><img src="/icons/back.gif" alt="[DIR]"></td>
	<td><a href="/gxfer3/">Parent Directory</a></td><td>&nbsp;</td>
	<td align="right">  - </td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td>
	<td><a href="Vervet_README.xls">Vervet_README.xls</a></td><td align="right">26-Jun-2012 14:06  </td>
	<td align="right"> 34K</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td>
	<td><a href="gerald_C0RCYACXX_5_AGATAG.bam">gerald_C0RCYACXX_5_AGATAG.bam</a></td>
	<td align="right">26-Jun-2012 09:38  </td><td align="right">3.3G</td></tr>
<tr><td valign="top"><img src="/icons/unknown.gif" alt="[   ]"></td>
	<td><a href="gerald_C0RCYACXX_5_AGATAG.bam.md5">gerald_C0RCYACXX_5_AGATAG.bam.md5</a></td>
	<td align="right">26-Jun-2012 14:06  </td><td align="right"> 64 </td></tr>
	...
	
				"""
				HTMLParser.__init__(self)
			
			def handle_starttag(self, tag, attrs):
				if tag=='a':
					self.startTagList.append(tag)
					for attr in attrs:
						if attr[0]=='href':
							self.hrefValueList.append(attr[1])
						
			def handle_endtag(self, tag):
				#self.endTagList.append(tag)
				pass
			def handle_data(self, data):
				"""
				if the last start tag is "a" (link) and the last href value is equal to this data
					then this data is one of the files to be downloaded.
				"""
				if len(self.startTagList)>0 and self.startTagList[-1]=='a' and self.hrefValueList[-1]==data:
					self.relativePathList.append(data)
				#print "Encountered some data  :", data
		parser = MyHTMLParser()
		parser.feed(url_content)
		
		sys.stderr.write(' found %s files/folders.\n'%(len(parser.relativePathList)))
		return parser.relativePathList
		
	
	def addJobs(self, inputURL=None, relativePathList =[], outputDir="", username=None, password=None, \
			transferOutput=True):
		"""
		2012.6.27
		"""
		
		sys.stderr.write("Adding wget jobs for %s input ... "%(len(relativePathList)))
		no_of_jobs= 0
		
		topOutputDir = outputDir
		topOutputDirJob = yh_pegasus.addMkDirJob(self, mkdir=self.mkdirWrap, outputDir=topOutputDir)
		no_of_jobs += 1
		returnData = PassingData()
		returnData.jobDataLs = []
		
		for relativePath in relativePathList:
			#2013.06.26 remove all "/" from  relativePath in case it's a folder
			relativePathNoFolder = relativePath.replace('/', '_')
			logFile = File('%s.log'%(relativePathNoFolder))
			wgetJob = self.addWgetJob(executable=self.wget, url=inputURL, relativePath=relativePath, \
						username=username, password=password,\
						targetFolder=outputDir, logFile=logFile, cut_dir_number=self.cut_dir_number, parentJobLs=[topOutputDirJob], extraDependentInputLs=[], \
						transferOutput=transferOutput, \
						extraArguments=None, job_max_memory=50)
			#include the tfam (outputList[1]) into the fileLs
			returnData.jobDataLs.append(PassingData(jobLs=[wgetJob], file=wgetJob.output, \
											fileLs=wgetJob.outputLs))
			no_of_jobs += 1
		sys.stderr.write("%s jobs.\n"%(no_of_jobs))
		
		return returnData	
	
	def registerCustomExecutables(self):
		"""
		2012.6.27
		"""
		namespace = self.namespace
		version = self.version
		operatingSystem = self.operatingSystem
		architecture = self.architecture
		clusters_size = self.clusters_size
		site_handler = self.site_handler
		
		executableList = []
		wget = Executable(namespace=namespace, name="wget", version=version, os=operatingSystem, \
							arch=architecture, installed=True)
		wget.addPFN(PFN("file://" + os.path.join(self.pathToWget), site_handler))
		executableList.append(wget)
		
		for executable in executableList:
			executable.addProfile(Profile(Namespace.PEGASUS, key="clusters.size", value="%s"%self.clusters_size))
			self.addExecutable(executable)
			setattr(self, executable.name, executable)
	
	def run(self):
		"""
		2012.6.27
		"""
		
		if self.debug:
			import pdb
			pdb.set_trace()
		
		self.registerJars()
		self.registerExecutables()
		self.registerCustomExecutables()
		
		relativePathList = self.getFilenamesToBeDownloaded(url=self.inputURL, username=self.username, password=self.password)
		#one file needs to registered so that replica catalog is not empty
		#but this file doesn't need to be actually used by any job.
		wgetFile = self.registerOneInputFile(inputFname=self.pathToWget)
		
		self.addJobs(inputURL=self.inputURL, relativePathList=relativePathList, outputDir=self.pegasusFolderName, username=self.username, \
					password=self.password, \
					transferOutput=True)
		outf = open(self.outputFname, 'w')
		self.writeXML(outf)

if __name__ == '__main__':
	main_class = WgetWorkflow
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	instance = main_class(**po.long_option2value)
	instance.run()
