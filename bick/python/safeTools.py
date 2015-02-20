import lsst.pipe.base  as pipeBase
import lsst.pex.config as pexConfig
import lsst.daf.persistence as dafPersist
import lsst.pipe.tasks.coaddBase as coaddBase


class ReadOnlyWrapper(object):
    """
    ReadOnlyWrapper

    This is intended to wrap a Butler, or butler-like
    object (such as a dataRef) to make it read-only.

    It can be used as follows:

       realButler = dafPersist.Butler()
       roButler = ReadOnlyWrapper(realButler)

       # this does nothing
       roButler.put("stuff", obj)

    Similarly:
    
       roDataRef = ReadOnlyWrapper(dataRef)

       # this also does nothing
       roDataRef.put("stuff", obj)
    
       # an attempt to get a real butler with getButler() also returns a ReadOnly object.
       roButler = roDataRef.getButler()
       roButler.put("stuff", obj)
    
    """
    
    # Methods to disable so they do nothing
    noops = ["put"]
    # Methods to permit, but to wrap the result in a ReadOnlyWrapper
    wrap  = ["getButler"]

    
    def __init__(self, thing):
        self._thing = thing
        
    def __getattr__(self, method):
        
        # if we're ignoring this method, return a null method which prints a warning
        if method in self.noops:
            def noop(*args, **kwargs):
                print "Ignoring %s() method. I'm a ReadOnly %s." % (method, str(self._thing))
                print "args: ", args
                print "kwargs: " , kwargs
            return noop
            
        if method in self.wrap:
            def wrap(*args, **kwargs):
                return ReadOnlyWrapper(getattr(self._thing, method)(*args, **kwargs))
            return wrap
            
        # otherwise, return the wrapped method
        return getattr(self._thing, method)

    def __getstate__(self):
        state = self.__dict__.copy()
        return state
    def __setstate__(self, dict):
        self.__dict__.update(dict)

        
class SafeButler(ReadOnlyWrapper):
    """Create a safe butler which is ReadOnlyWrapper'd """
    def __init__(self, *args, **kwargs):
        super(SafeButler, self).__init__(dafPersist.Butler(*args, **kwargs))



        
class SafeTaskRunner(pipeBase.TaskRunner):
    """ A Special TaskRunner which makes sure that the dataRefs which are passed
    to the task's run() method are rendered read-only.

    All that's needed is to intercept the dataRef in getTargetList() and make it read-only.
    This is before anyone will get to use it.
    """
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(ReadOnlyWrapper(ref), kwargs) for ref in parsedCmd.id.refList]

class SafeConfig(pexConfig.Config):
    coaddName = pexConfig.ChoiceField(
        dtype   = str,
        doc     = "Type of coadd to use",
        default = "deep",
        allowed = {"deep"   : "deepCoadd"}
        )

butlerTarget = "raw"
dataIdContainer  = {
    "raw":              pipeBase.DataIdContainer,
    "src":              pipeBase.DataIdContainer,
    "deepCoadd":        coaddBase.ExistingCoaddDataIdContainer,
    "deepCoadd_skyMap": coaddBase.ExistingCoaddDataIdContainer,
    "wcs":              coaddBase.ExistingCoaddDataIdContainer,
}
    
class SafeTask(pipeBase.CmdLineTask):
    """A light-weight utility task with a read-only butler that hides
    some of the boiler-plate of task writing.
    
    Derived objects are intended to be used for reading and working with data products,
    rather than acting as components of the pipeline.  The butler/dataRefs are wrapped to
    be read-only to prevent accidentally overwriting data on disk.

    The user can derive from this class and need only specify scatter and gather
    (both are optional).
    """

    _DefaultName = "Safe"
    ConfigClass  = SafeConfig
    RunnerClass  = SafeTaskRunner

    @classmethod
    def _makeArgumentParser(cls):
        parser = pipeBase.ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", butlerTarget, help="Data ID, e.g. --id tract=1234 patch=2,2",
                               ContainerClass=dataIdContainer[butlerTarget])
        return parser
    
    @classmethod
    def parseAndRun(cls, *args, **kwargs):
        setattr(cls, "_getConfigName",       SafeTask._getConfigName)
        setattr(cls, "_getEupsVersionsName", SafeTask._getEupsVersionsName)
        setattr(cls, "_getMetadataName",     SafeTask._getMetadataName)
        kwargs['doReturnResults'] = True
        results = super(SafeTask,cls).parseAndRun(*args, **kwargs)
        all = []
        for result in results.resultList:
            all.append(result.result)
        cls().gather(all)

    def run(self, *args, **kwargs):
        pass
    def gather(self, results):
        pass
        
    # we need to kill each of these methods so the users won't
    # mess with persisted configs and version info
    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None


        
