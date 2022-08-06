## Specifics for NAF-CMS
Problem was that I could not get the submission to work when sourcing crab3 via `source crab.sh`, despite having activated grid-UI, valid proxy, in cmsenv (all using `bash` shell). Setup was sourcing grid-commands, cmsset and crab already in `.bashrc`. Testing worked perfectly fine locally.

<details>
  <summary>Expand for error message using current instructions for CRAB-submission with PFNano on DESY NAF</summary>
    
```
==> /DoubleMuon/Run2017B-09Aug2019_UL2017-v1/MINIAOD
Process Process-5:
Traceback (most recent call last):
  File "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/python/2.7.15-pafccj/lib/python2.7/multiprocessing/process.py", line 267, in _bootstrap
    self.run()
  File "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/python/2.7.15-pafccj/lib/python2.7/multiprocessing/process.py", line 114, in run
    self._target(*self._args, **self._kwargs)
  File "crabby.py", line 16, in submit
    crabCommand('submit', config = config)
  File "/cvmfs/cms.cern.ch/share/cms/crab-prod/v3.220323.01/lib/CRABAPI/RawCommand.py", line 28, in crabCommand
    return execRaw(command, arguments)
  File "/cvmfs/cms.cern.ch/share/cms/crab-prod/v3.220323.01/lib/CRABAPI/RawCommand.py", line 47, in execRaw
    cmdobj = getattr(mod, command)(logger, args)
  File "/cvmfs/cms.cern.ch/share/cms/crab-prod/v3.220323.01/lib/CRABClient/Commands/submit.py", line 36, in __init__
    SubCommand.__init__(self, logger, cmdargs, disable_interspersed_args=True)
  File "/cvmfs/cms.cern.ch/share/cms/crab-prod/v3.220323.01/lib/CRABClient/Commands/SubCommand.py", line 375, in __init__
    self.handleMyProxy()  
  File "/cvmfs/cms.cern.ch/share/cms/crab-prod/v3.220323.01/lib/CRABClient/Commands/SubCommand.py", line 492, in handleMyProxy
    raise ProxyCreationException("Problems delegating My-proxy.\n%s" % msg1)
ProxyCreationException: Problems delegating My-proxy.
Error trying to create credential:
 Problems delegating My-proxy. Error executing export GT_PROXY_MODE=rfc ; myproxy-init -d -n -s myproxy.cern.ch -C ~/.globus/usercert.pem -y ~/.globus/userkey.pem -x -R '/DC=ch/DC=cern/OU=computers/CN=crab-(preprod|prod|dev
)-tw(01|02|03).cern.ch|/DC=ch/DC=cern/OU=computers/CN=stefanov(m|m2).cern.ch|/DC=ch/DC=cern/OU=computers/CN=dciangot-tw.cern.ch' -x -Z '/DC=ch/DC=cern/OU=computers/CN=crab-(preprod|prod|dev)-tw(01|02|03).cern.ch|/DC=ch/DC=c
ern/OU=computers/CN=stefanov(m|m2).cern.ch|/DC=ch/DC=cern/OU=computers/CN=dciangot-tw.cern.ch' -l 887641ffe656e41de4f4f2e594f61348ad064be7 -t 168 -c 720:00:

myproxy-init: error while loading shared libraries: libmyproxy.so.6: cannot open shared object file: No such file or directory

```
</details>
<br>

Somehow fixed when switching to different ui-profile & crab-setup (not sure if it's only a coincidence because something else was wrong I am not aware of):
```
source /cvmfs/grid.cern.ch/centos7-umd4-ui-4_200423/etc/profile.d/setup-c7-ui-example.sh
source /cvmfs/cms.cern.ch/common/crab-setup.sh prod
source /cvmfs/cms.cern.ch/cmsset_default.sh
< navigate to CMSSW_X_Y_Z/src >
cmsenv
```
Used both `grid-proxy-init` as well as `voms-proxy-init -vms cms`, yet was asked upon submission to type in password for certificate.