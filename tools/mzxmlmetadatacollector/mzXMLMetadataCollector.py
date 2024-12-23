# -*- coding: utf-8 -*-
"""
Created on Dec 04 2023
@author: quruin

Last update on Nov 14 2024

@author: quruin
"""
# Packages
import os
import re
import sys

import numpy as np


def get_info(file, format1):
    # mzml files #
    if (format1.casefold() == 'mzml'):
        # Saving all info of first file until first scan
        f = open(file, 'r+', encoding="utf-8")
        with f:
            text1 = ''
            while ('<binary>' not in text1):
                text1 = text1 + f.readline()
        f.close()
        # Getting version format
        try:
            version1 = text1.split('<mzML')[1] \
                .split('version="')[1].split('"')[0]
        except Exception:
            version1 = "Not found"
            raise
        # Getting ko size
        try:
            taillek1 = np.round(os.stat(file).st_size/1024, 2)
        except Exception:
            taillek1 = "Calculation failed"
            raise
        # Getting Mo size
        try:
            taillem1 = np.round(os.stat(file).st_size/1024**2, 2)
        except Exception:
            taillem1 = "Calculation failed"
            raise
        # Getting Go size
        try:
            tailleg1 = np.round(os.stat(file).st_size/1024**3, 2)
        except Exception:
            tailleg1 = "Calculation failed"
            raise
        # Getting MS level
        try:
            if ('<cvParam cvRef="MS" accession="MS:1000580" name="' in text1):
                mslevel1 = text1.split('<cvParam cvRef="MS" \
                    accession="MS:1000580" name="')[1].split('"')[0]
            else:
                if ('<cvParam cvRef="MS" accession="MS:1000580" \
                        value="" name="' in text1):
                    mslevel1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1000580" value="" name="')[1]\
                            .split('"')[0]
                else:
                    if ('accession="MS:1000579" cvRef="MS" \
                            name="' in text1):
                        mslevel1 = text1.split('accession="MS:1000579" \
                            cvRef="MS" name="')[1].split('"')[0]
                    else:
                        if ('cvRef="MS" accession="MS:1000579" \
                                name="' in text1):
                            mslevel1 = text1.split('cvRef="MS" \
                                accession="MS:1000579" name="')[1]\
                                    .split('"')[0]
                        else:
                            if (' <cvParam cvRef="MS" \
                                    accession="MS:1000580" name="' in text1):
                                mslevel1 = text1.split(' <cvParam cvRef="MS" \
                                    accession="MS:1000580" name="')[1]\
                                        .split('"')[0]
                            else:
                                mslevel1 = "Not found"
        except Exception:
            mslevel1 = 'Error'
            raise
        # Getting spectrum type
        try:
            if ('<cvParam cvRef="MS" accession="MS:1000127" \
                    name="' in text1):
                spectrum1 = text1.split('<cvParam cvRef="MS"  \
                    accession="MS:1000127" name="')[1].split('"')[0]
            else:
                if ('<cvParam cvRef="MS" accession="MS:1000127" \
                        value="" name="' in text1):
                    spectrum1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1000127" value="" name="')[1]\
                            .split('"')[0]
                else:
                    if ('<cvParam accession="MS:1000127" \
                            cvRef="MS" name="' in text1):
                        spectrum1 = text1.split('<cvParam \
                            accession="MS:1000127" \
                                cvRef="MS" name="')[1].split('"')[0]
                    else:
                        spectrum1 = "Not found"
        except Exception:
            spectrum1 = "Not found"
            raise
        # Getting sourc file
        try:
            source1 = text1.split('<sourceFile ')[1].split('">')[0]\
                .split('location="')[1].split('"')[0]
        except Exception:
            try:
                source1 = text1.split('<sourceFile id="WIFF" name="')[1]\
                    .split('"')[0]
            except Exception:
                try:
                    source1 = text1.split('<sourceFile id=" ')[1].split('"')[0]
                except Exception:
                    source1 = "Not found"
                    raise
        # Getting acquisition date
        try:
            date1 = text1.split('startTimeStamp="')[1].split('"')[0]\
                .split('T')[0] + " " + text1.split('startTimeStamp="')[1]\
                .split('"')[0].split('T')[1]
        except Exception:
            date1 = "Not found"
            raise
        # Getting software used
        try:
            softwaresList1 = ''
            subtext1 = text1.split('softwareList count')[1] \
                .split('</softwareList>')[0]
            for i in range(subtext1.count("<software")):
                if i != 0:
                    softwaresList1 = softwaresList1 + ' + '
                softwaresversions1 = subtext1.split("<software")[i+1]\
                    .split('version="')[1].split('"')[0]
                if ('cvRef="MS" name="' in subtext1):
                    softwares1 = subtext1.split('cvRef="MS" name="')[i+1]\
                        .split('"')[0]
                else:
                    softwares1 = subtext1.split('<software id="')[i+1]\
                        .split('"')[0]
                softwaresList1 = softwaresList1 + \
                    softwares1 + ' ' + softwaresversions1
        except Exception:
            softwaresList1 = "Not found"
            raise
        # Getting processing methods
        try:
            ProcessList1 = ''
            subtext1 = text1.split('<dataProcessingList count')[1] \
                .split('</dataProcessingList>')[0]
            for i in range(subtext1.count("<processingMethod ")):
                if i != 0:
                    ProcessList1 = ProcessList1 + ' + '
                methods1 = subtext1.split("<processingMethod ")[i+1] \
                    .split('name="')[1].split('"')[0]
                softwares1 = subtext1.split('softwareRef="')[1] \
                    .split('"')[0]
                ProcessList1 = ProcessList1 + softwares1 + ' ' + methods1
        except Exception:
            ProcessList1 = "Not found"
            raise
        # Getting machine model
        try:
            if ('<cvParam cvRef="MS" accession="MS:1000703" name="' in text1):
                modele1 = text1.split('<cvParam cvRef="MS" \
                    accession="MS:1000703" name="')[1].split('"')[0]
            else:
                if ('accession="MS:1001547"' in text1):
                    modele1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1001547" name="')[1].split('"')[0]
                else:
                    if ('<cvParam cvRef="MS" \
                            accession="MS:1003094" name="' in text1):
                        modele1 = text1.split('<cvParam cvRef="MS" \
                            accession="MS:1003094" name="')[1].split('"')[0]
                    else:
                        if ('accession="MS:1003123" name="' in text1):
                            modele1 = text1.split('<cvParam cvRef="MS" \
                                accession="MS:1003123" name="')[1] \
                                    .split('"')[0]
                        else:
                            if ('<cvParam cvRef="MS" \
                                    accession="MS:1000495"' in text1):
                                modele1 = text1.split('<cvParam cvRef="MS" \
                                    accession="MS:1000495" name="')[1]\
                                        .split('value="')[1].split('"')[0]
                            else:
                                if ('<cvParam accession="MS:1000703" \
                                        cvRef="MS" name="' in text1):
                                    modele1 = text1.split('<cvParam \
                                        accession="MS:1000703" \
                                            cvRef="MS" name="')[1] \
                                                .split('"')[0]
                                else:
                                    if ('<cvParam cvRef="MS" accession= \
                                        "MS:1000483" \
                                            value="" name="' in text1):
                                        modele1 = text1.split('<cvParam \
                                            cvRef="MS" accession="MS:1000483" \
                                                value="" name="')[1] \
                                                    .split('"')[0]
                                    else:
                                        if ('<cvParam accession="MS:1000703" \
                                                cvRef="MS" name="' in text1):
                                            modele1 = text1.split('<cvParam \
                                                accession="MS:1000703" \
                                                    cvRef="MS" \
                                                        name="')[1] \
                                                            .split('"')[0]
                                        else:
                                            modele1 = 'Not available'
            if (modele1 == ''):
                modele1 = 'Not available'
        except Exception:
            modele1 = "Not found"
            raise
        # Getting points and scans number
        try:
            nbscans1 = text1.split('<spectrumList count="')[1].split('"')[0]
        except Exception:
            nbscans1 = "Not found"
            raise
        try:
            nbpoints1 = 0
            with open(file, 'r') as f:
                line = 'start'
                while (line != ''):
                    line = f.readline()
                    if ('defaultArrayLength=' in line):
                        nbpoints1 = nbpoints1 + int(line
                        .split('defaultArrayLength="')[1]
                        .split('"')[0])
            f.close()
        except Exception:
            nbpoints1 = "Calculation failed"
            raise
        # Encoding
        encodage1 = re.search('name="(.*)-bit float"', text1) \
            .group(0).replace('name="', "").replace('"', '')
        return format1, version1, taillek1, taillem1, tailleg1, mslevel1, \
            spectrum1, source1, date1, softwaresList1, ProcessList1, \
                modele1, nbscans1, nbpoints1, encodage1
    # mzXml files #
    if (format1.casefold() == 'mzxml'):
        # Saving all info of first file until first scan
        f = open(file, 'r+', encoding="utf-8")
        with f:
            text1 = ''
            while ('</scan>' not in text1):
                text1 = text1 + f.readline()
        f.close()
        # Getting format version
        try:
            version1 = text1.split('<mzXML')[1].split('mzXML_')[1]\
                .split('"')[0].split(' ')[0]
        except Exception:
            version1 = "Not found"
            raise
        # Getting ko size
        try:
            taillek1 = np.round(os.stat(file).st_size/1024, 2)
        except Exception:
            taillek1 = "Calculation failed"
            raise
        # Getting Mo size
        try:
            taillem1 = np.round(os.stat(file).st_size/1024**2, 2)
        except Exception:
            taillem1 = "Calculation failed"
            raise
        # Getting Go size
        try:
            tailleg1 = np.round(os.stat(file).st_size/1024**3, 2)
        except Exception:
            tailleg1 = "Calculation failed"
            raise
        # Getting MS level
        mslevel1 = "Not available in mzXML"
        # Getting spectrum type
        try:
            spectrum1 = text1.split('centroided="')[1].split('"')[0]
            if spectrum1 == "0":
                spectrum1 = "profile"
            else:
                if spectrum1 == "1":
                    spectrum1 = "centroid"
                else:
                    spectrum1 = "Not found"
        except Exception:
            spectrum1 = "Not found"
            raise
        # Getting source file
        try:
            source1 = text1.split('parentFile fileName="')[1].split('"')[0]
        except Exception:
            source1 = "Not found"
            raise
        # Getting acquisition date
        date1 = "Not available in mzXML"
        # Getting used softwares
        try:
            softwaresList1 = ''
            for i in range(text1.count("<software")):
                if i != 0:
                    softwaresList1 = softwaresList1 + ' + '
                softwaresList1 = softwaresList1 + \
                    text1.split("<software")[i+1].split('type="')[1]\
                    .split('"')[0] + ': ' +\
                        text1.split("<software")[i+1].split('name="')[1]\
                            .split('"')[0] + ' ' +\
                                text1.split("<software")[i+1]\
                                    .split('version="')[1]\
                                        .split('"')[0]
        except Exception:
            softwaresList1 = "Not found"
            raise
        # Processing
        ProcessList1 = "Not available in mzXML"
        # Getting Machine model
        try:
            if ('<msModel category="msModel" value="' in text1):
                modele1 = text1.split('<msModel category="msModel" \
                    value="')[1].split('"')[0]
            else:
                if ('accession="MS:1001547"' in text1):
                    modele1 = text1.split('<cvParam cvRef="MS" \
                        accession="MS:1001547" name="')[1].split('"')[0]
                else:
                    if ('accession="MS:1003123"' in text1):
                        modele1 = text1.split('<cvParam cvRef="MS" \
                            accession="MS:1003123" name="')[1].split('"')[0]
                    else:
                        if ('accession="MS:1000495"' in text1):
                            modele1 = text1.split('<cvParam cvRef="MS" \
                                accession="MS:1000495" name="')[1]\
                                    .split('"')[0]
                        else:
                            modele1 = "Not found"
        except Exception:
            modele1 = "Not found"
            raise
        # Getting scans and points number
        try:
            nbscans1 = text1.split('msRun scanCount="')[1].split('"')[0]
        except Exception:
            nbscans1 = "Not found"
            raise
        try:
            nbpoints1 = 0
            with open(file, 'r') as f:
                line = 'start'
                while (line != ''):
                    line = f.readline()
                    if ('peaksCount="' in line):
                        nbpoints1 = nbpoints1 + \
                            int(line.split('peaksCount="')[1].split('"')[0])
        except Exception:
            nbpoints1 = "Calculation failed"
            raise
        # Encoding
        try:
            encodage1 = text1.split('precision="')[1].split('"')[0] + "-bit"
        except Exception:
            encodage1 = 'Not available'
            raise
        return format1, version1, taillek1, taillem1, tailleg1, \
            mslevel1, spectrum1, source1, date1, softwaresList1, \
            ProcessList1, \
            modele1, nbscans1, nbpoints1, encodage1
    # Other types of files
    ######################################################
    if ((format1.casefold() != 'mzml') & (format1.casefold() != 'mzxml')):
        return ['' for i in range(0, 15)]


def update_entries(infile, singleFile):
    if (singleFile == "SINGLE"):
        info = ''
        info += 'Name\tformat\tversion\tsize(ko)\tsize(Mo)\tsize(Go)\
        \tMSlevel\tSpectrum type\tSource file\tAcquisition \
        date\tSoftware(s) used\tProcessing method(s)\tMachine\t \
        Number of scans\tNumber of points\tEncoding\n'
        x = get_info(infile, sys.argv[4].split('.')[-1])
        info += sys.argv[4] + '\t'
        for i in range(len(x)):
            info += str(x[i]) + '\t'
        info += '\n'
        return info
    else:
        info = ''
        info += 'Name\tformat\tversion\tsize(ko)\tsize(Mo)\tsize(Go)\t \
        MSlevel\tSpectrum type\tSource file\tAcquisition \
        date\tSoftware(s) used\tProcessing method(s)\tMachine\t \
        Number of scans\tNumber of points\tEncoding\n'
        ii = 4
        for f in infile.split(','):
            x = get_info(f, sys.argv[ii].split('.')[-1])
            info += sys.argv[ii] + '\t'
            for i in range(len(x)):
                info += str(x[i]) + '\t'
            info += '\n'
            ii += 1
        return info


print([i for i in sys.argv])
outfile = sys.argv[1]
singleFile = sys.argv[2]
infile = sys.argv[3]
res = open(outfile, "w")
res.write(update_entries(infile, singleFile))
res.close()
